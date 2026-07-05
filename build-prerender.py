#!/usr/bin/env python3
"""
build-prerender.py
讀取 papers-data.js、projects-data.js 和 shared.js，產生預渲染 HTML，
並自動嵌入 index.html 的對應容器中。
同時產生 JSON-LD 結構化資料供搜尋引擎使用。

預渲染範圍：
  - 著作（papers-data.js → #papers-container）
  - 計畫（projects-data.js → #projects-container）
  - 自傳（shared.js bioZh/bioEn → #bio-content）
  - 榮譽（shared.js honorsZh/honorsEn → #honors-list）
  - 傑出校友推薦文（shared.js alumniArticleZh/alumniArticleEn → #alumni-article）
  - 杏壇芬芳推薦文（shared.js honorsArticleZh/honorsArticleEn → #honors-article）
  - 統計卡（stats-data.js statsZh/statsEn → #stats-container）
  - 導覽列（shared.js navItemsZh/navItemsEn → #main-nav，無 JS 後備）

用法：python3 build-prerender.py [--bump-dates]
  --bump-dates  同步更新五處日期戳記（DCTERMS.modified、sitemap 首頁
                lastmod、humans.txt、llms-full.txt、README 最後更新）；
                內容有實質變動時使用。
任何解析失敗、容器遺失或數量低於門檻都會以非零狀態退出且不寫入 index.html。
"""

import json
import os
import re
import sys

BASE = os.path.dirname(os.path.abspath(__file__))

# ─── Fail-fast error collector ───
# 任何解析失敗或容器遺失都記錄於此；main() 會在寫入 index.html 之前
# 檢查並以非零狀態退出，避免以空內容覆蓋正確的預渲染結果。
BUILD_ERRORS = []


def build_error(msg):
    print(f"❌ {msg}")
    BUILD_ERRORS.append(msg)


# ─────────────────────────────────────────────
# Parse JavaScript data using a bracket-based approach
# ─────────────────────────────────────────────

def extract_js_string(text, start, quote="'"):
    """Extract a JS quoted string starting at position `start` (must be the opening
    quote; single or double). Returns (string_value, end_position_after_closing_quote)."""
    assert text[start] == quote, f"Expected {quote} at pos {start}, got {text[start]!r}"
    i = start + 1
    chars = []
    while i < len(text):
        c = text[i]
        if c == '\\' and i + 1 < len(text):
            nxt = text[i + 1]
            if nxt == 'u' and i + 5 < len(text):
                code = int(text[i + 2:i + 6], 16)
                i += 6
                if 0xD800 <= code <= 0xDBFF and text[i:i + 2] == '\\u':
                    low = int(text[i + 2:i + 6], 16)
                    code = 0x10000 + ((code - 0xD800) << 10) + (low - 0xDC00)
                    i += 6
                chars.append(chr(code))
            elif nxt == 'n':
                chars.append('\n')
                i += 2
            elif nxt == 't':
                chars.append('\t')
                i += 2
            else:
                chars.append(nxt)
                i += 2
        elif c == quote:
            return ''.join(chars), i + 1
        else:
            chars.append(c)
            i += 1
    raise ValueError("Unterminated string")


def find_bracket_end(text, start, open_b='[', close_b=']'):
    """Find the matching closing bracket, skipping single- and double-quoted strings."""
    assert text[start] == open_b
    depth = 1
    i = start + 1
    while i < len(text) and depth > 0:
        c = text[i]
        if c == "'" or c == '"':
            _, i = extract_js_string(text, i, c)
        elif c == open_b:
            depth += 1
            i += 1
        elif c == close_b:
            depth -= 1
            i += 1
        else:
            i += 1
    return i  # position after closing bracket


def parse_papers_data(js_text):
    """Parse papersData array from papers-data.js"""
    categories = []

    # Find the top-level array start
    m = re.search(r'const\s+papersData\s*=\s*\[', js_text)
    if not m:
        build_error("Could not find papersData in papers-data.js")
        return categories

    arr_start = m.end() - 1  # position of '['
    arr_end = find_bracket_end(js_text, arr_start)
    arr_content = js_text[arr_start + 1:arr_end - 1]

    # Find each category object (top-level { ... } inside the array)
    i = 0
    while i < len(arr_content):
        # Find next {
        idx = arr_content.find('{', i)
        if idx == -1:
            break
        obj_end = find_bracket_end(arr_content, idx, '{', '}')
        obj_text = arr_content[idx:obj_end]
        i = obj_end

        cat = {'titleZh': '', 'titleEn': '', 'noteZh': None, 'noteEn': None, 'periods': []}

        # Extract titleZh
        tm = re.search(r"titleZh:\s*'", obj_text)
        if tm:
            cat['titleZh'], _ = extract_js_string(obj_text, tm.end() - 1)

        # Extract titleEn
        tm = re.search(r"titleEn:\s*'", obj_text)
        if tm:
            cat['titleEn'], _ = extract_js_string(obj_text, tm.end() - 1)

        # Extract noteZh (optional)
        tm = re.search(r"noteZh:\s*'", obj_text)
        if tm:
            cat['noteZh'], _ = extract_js_string(obj_text, tm.end() - 1)

        # Extract noteEn (optional)
        tm = re.search(r"noteEn:\s*'", obj_text)
        if tm:
            cat['noteEn'], _ = extract_js_string(obj_text, tm.end() - 1)

        # Find periods array
        pm = re.search(r'periods:\s*\[', obj_text)
        if pm:
            periods_start = pm.end() - 1
            periods_end = find_bracket_end(obj_text, periods_start)
            periods_content = obj_text[periods_start + 1:periods_end - 1]

            # Parse each period object
            j = 0
            while j < len(periods_content):
                pidx = periods_content.find('{', j)
                if pidx == -1:
                    break
                pobj_end = find_bracket_end(periods_content, pidx, '{', '}')
                pobj_text = periods_content[pidx:pobj_end]
                j = pobj_end

                period = {'titleZh': None, 'titleEn': None, 'items': []}

                # Extract titleZh
                ptm = re.search(r"titleZh:\s*'", pobj_text)
                if ptm:
                    period['titleZh'], _ = extract_js_string(pobj_text, ptm.end() - 1)
                # Check for null
                ptm_null = re.search(r"titleZh:\s*null", pobj_text)

                # Extract titleEn
                ptm = re.search(r"titleEn:\s*'", pobj_text)
                if ptm:
                    period['titleEn'], _ = extract_js_string(pobj_text, ptm.end() - 1)

                # Find items array
                im = re.search(r'items:\s*\[', pobj_text)
                if im:
                    items_start = im.end() - 1
                    items_end = find_bracket_end(pobj_text, items_start)
                    items_content = pobj_text[items_start + 1:items_end - 1]

                    # Extract each string item
                    k = 0
                    while k < len(items_content):
                        # Find next quote
                        qidx = items_content.find("'", k)
                        if qidx == -1:
                            break
                        item_str, k = extract_js_string(items_content, qidx)
                        period['items'].append(item_str)

                cat['periods'].append(period)

        categories.append(cat)

    return categories


def parse_projects_data(js_text):
    """Parse projectsData array from projects-data.js"""
    roles = []

    m = re.search(r'const\s+projectsData\s*=\s*\[', js_text)
    if not m:
        build_error("Could not find projectsData in projects-data.js")
        return roles

    arr_start = m.end() - 1
    arr_end = find_bracket_end(js_text, arr_start)
    arr_content = js_text[arr_start + 1:arr_end - 1]

    i = 0
    while i < len(arr_content):
        idx = arr_content.find('{', i)
        if idx == -1:
            break
        # This is a role object: { titleZh, titleEn, items: [...] }
        obj_end = find_bracket_end(arr_content, idx, '{', '}')
        obj_text = arr_content[idx:obj_end]
        i = obj_end

        role = {'titleZh': '', 'titleEn': '', 'items': []}

        tm = re.search(r"titleZh:\s*'", obj_text)
        if tm:
            role['titleZh'], _ = extract_js_string(obj_text, tm.end() - 1)

        tm = re.search(r"titleEn:\s*'", obj_text)
        if tm:
            role['titleEn'], _ = extract_js_string(obj_text, tm.end() - 1)

        # Find items array
        im = re.search(r'items:\s*\[', obj_text)
        if im:
            items_start = im.end() - 1
            items_end = find_bracket_end(obj_text, items_start)
            items_content = obj_text[items_start + 1:items_end - 1]

            # Parse each project item object
            j = 0
            while j < len(items_content):
                pidx = items_content.find('{', j)
                if pidx == -1:
                    break
                pobj_end = find_bracket_end(items_content, pidx, '{', '}')
                pobj_text = items_content[pidx:pobj_end]
                j = pobj_end

                proj = {'period': '', 'source': '', 'grant': '', 'internal': None, 'amount': None, 'name': ''}

                for key in ['period', 'source', 'grant', 'internal', 'amount', 'name']:
                    km = re.search(rf"{key}:\s*'", pobj_text)
                    if km:
                        proj[key], _ = extract_js_string(pobj_text, km.end() - 1)
                    else:
                        km_null = re.search(rf"{key}:\s*null", pobj_text)
                        if km_null:
                            proj[key] = None

                role['items'].append(proj)

        roles.append(role)

    return roles


# ─── HTML escaping ───

def esc(s):
    if s is None:
        return ''
    return s.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;').replace('"', '&quot;')


# ─── Replicate formatPaperText from shared.js ───

def format_paper_text(text):
    html = text.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')

    html = re.sub(r'Huang C\.Y\.(\*?)', r'<strong>Huang C.Y.\1</strong>', html)
    html = re.sub(r'黃崇源', '<strong>黃崇源</strong>', html)

    html = re.sub(r'(https?://[^\s,)]+[^\s,).])', r'<a href="\1" target="_blank" rel="noopener noreferrer">\1</a>', html)
    html = re.sub(r'doi:(10\.\d{4,}/[^\s,)]+[^\s,).])', r'<a href="https://doi.org/\1" target="_blank" rel="noopener noreferrer">doi:\1</a>', html)
    html = re.sub(r'\[PDF:([^\]]+)\]', r'<a href="\1" target="_blank" rel="noopener noreferrer">[PDF]</a>', html)

    html = re.sub(r'\(([^()]*(?:SCI|SSCI|SCIE|EI)[^()]*)\)\s*\.?\s*$', r'<span class="paper-index">(\1)</span>', html)

    return html


# ─── Idempotent container replacement ───

def replace_container_inner(html, open_tag, new_inner, tag_name):
    """Replace the inner content of a container, whether it is empty or already
    filled by a previous run, by depth-matching nested tags of the same name."""
    start = html.find(open_tag)
    if start < 0:
        return html, False
    inner_start = start + len(open_tag)
    open_re = re.compile('<' + tag_name + r'[\s>]')
    close_tag = '</' + tag_name + '>'
    depth, i = 1, inner_start
    while True:
        m = open_re.search(html, i)
        c = html.find(close_tag, i)
        if c < 0:
            return html, False
        if m and m.start() < c:
            depth += 1
            i = m.end()
        else:
            depth -= 1
            if depth == 0:
                return html[:inner_start] + new_inner + html[c:], True
            i = c + len(close_tag)


# ─── Generate Papers HTML ───

def generate_papers_html(papers_data):
    lines = []
    for category in papers_data:
        lines.append('<div class="papers-category">')
        lines.append(f'  <h3 class="papers-category-title">{esc(category["titleZh"])}</h3>')

        if category.get('noteZh'):
            lines.append(f'  <p class="papers-note">{esc(category["noteZh"])}</p>')

        global_index = 1

        for period in category['periods']:
            if period.get('titleZh'):
                lines.append(f'  <h4 class="papers-period-title">{esc(period["titleZh"])}</h4>')

            lines.append(f'  <ol class="papers-list" start="{global_index}">')

            for text in period['items']:
                lines.append(f'    <li class="paper-item">{format_paper_text(text)}</li>')
                global_index += 1

            lines.append('  </ol>')

        lines.append('</div>')

    return '\n'.join(lines)


# ─── Generate Projects HTML ───

def generate_projects_html(projects_data):
    lines = []

    for role in projects_data:
        lines.append('<div class="projects-role">')
        lines.append(f'  <h3 class="projects-role-title">{esc(role["titleZh"])}</h3>')
        lines.append('  <ol class="projects-list">')

        for proj in role['items']:
            parts = [esc(proj['source']), '期間：' + esc(proj['period'])]
            if proj.get('amount'):
                parts.append(esc(proj['amount']))

            grant_str = esc(proj['grant'])
            if proj.get('internal'):
                grant_str += ' / ' + esc(proj['internal'])

            lines.append('    <li class="project-item">')
            lines.append(f'      <span class="project-name">{esc(proj["name"])}</span>')
            lines.append(f'      <span class="project-meta">{" &middot; ".join(parts)}'
                         f'<br><span class="project-grant">{grant_str}</span></span>')
            lines.append('    </li>')

        lines.append('  </ol>')
        lines.append('</div>')

    return '\n'.join(lines)


# ─── Generate JSON-LD for publications ───

def generate_jsonld(papers_data, projects_data):
    articles = []
    skipped = []

    # 以分類標題選取期刊分類（不可依賴陣列位置——重排資料檔不應造成錯置）
    journal_cat = next((c for c in papers_data if 'SCI' in c['titleZh']), None)
    if papers_data and journal_cat is None:
        build_error("JSON-LD: no papers category with 'SCI' in titleZh found")

    if journal_cat:
        for period in journal_cat['periods']:
            for text in period['items']:
                year_match = re.search(r'\((\d{4})\)', text)
                year = year_match.group(1) if year_match else None

                title_match = re.search(
                    r'\(\d{4}\)\s+(.+?)\.\s+(?:IEEE|PLoS|Physica|Scientific|Computational|Journal|Simulation|'
                    r'International|Applied|BioMed|Artificial|Discrete|WSEAS|Mathematical|The\s+Journal)',
                    text
                )
                title = title_match.group(1) if title_match else None

                doi_match = re.search(r'doi:(10\.\d{4,}/[^\s,)]+[^\s,).])', text)
                doi = doi_match.group(1) if doi_match else None

                url_match = re.search(r'(https?://[^\s,)]+[^\s,).])', text)
                url = url_match.group(1) if url_match else (f'https://doi.org/{doi}' if doi else None)

                if title and year:
                    article = {
                        "@type": "ScholarlyArticle",
                        "name": title,
                        "datePublished": year,
                        "author": {"@type": "Person", "@id": "https://canslab1.github.io/#person"}
                    }
                    if doi:
                        article["sameAs"] = f"https://doi.org/{doi}"
                    if url:
                        article["url"] = url
                    articles.append(article)
                else:
                    skipped.append(text[:120])

        # 漏抓偵測：期刊白名單或年份 regex 沒抓到的論文會從結構化資料
        # 靜默消失（例如未來發表於 Nature/Frontiers 等不在白名單的期刊）
        cat_total = sum(len(p['items']) for p in journal_cat['periods'])
        if len(articles) != cat_total:
            build_error(f"JSON-LD dropped {cat_total - len(articles)} of {cat_total} journal papers "
                        f"(title/year regex miss — extend the journal-name whitelist at generate_jsonld)")
            for s in skipped:
                print(f"    · skipped: {s}")

    research_projects = []
    # 以角色標題選取主持人分類（不可依賴陣列位置；須精確比對，
    # 因「共同主持人」亦包含「主持人」子字串）
    pi_cat = next((r for r in projects_data if r['titleZh'] == '主持人'), None)
    if projects_data and pi_cat is None:
        build_error("JSON-LD: no projects role with titleZh == '主持人' found")
    if pi_cat:
        for proj in pi_cat['items']:
            rp = {
                "@type": "ResearchProject",
                "name": proj['name'],
                "funder": proj['source'],
                "identifier": proj['grant'],
                "member": {"@type": "Person", "@id": "https://canslab1.github.io/#person"}
            }
            if proj.get('amount'):
                rp["funding"] = proj['amount']
            research_projects.append(rp)

    items = []
    pos = 1
    for a in articles:
        items.append({"@type": "ListItem", "position": pos, "item": a})
        pos += 1
    for p in research_projects:
        items.append({"@type": "ListItem", "position": pos, "item": p})
        pos += 1

    ld = {
        "@context": "https://schema.org",
        "@type": "ItemList",
        "name": "Publications and Research Projects of Prof. Chung-Yuan Huang",
        "description": "Complete list of scholarly publications and funded research projects by Prof. Chung-Yuan Huang (黃崇源教授), Chang Gung University.",
        "numberOfItems": len(articles) + len(research_projects),
        "itemListElement": items
    }

    return ld


# ─── Parse simple string arrays from shared.js ───

def parse_string_array(js_text, var_name):
    """Parse a const array of strings from shared.js, e.g. const bioZh = ['...', '...']
    Returns a list of strings."""
    m = re.search(rf'const\s+{var_name}\s*=\s*\[', js_text)
    if not m:
        build_error(f"Could not find {var_name} in shared.js")
        return []

    arr_start = m.end() - 1  # position of '['
    arr_end = find_bracket_end(js_text, arr_start)
    arr_content = js_text[arr_start + 1:arr_end - 1]

    items = []
    i = 0
    while i < len(arr_content):
        qidx = arr_content.find("'", i)
        if qidx == -1:
            break
        item_str, i = extract_js_string(arr_content, qidx)
        items.append(item_str)

    return items


# ─── Generate Bio HTML ───

def generate_bio_html(bio_zh, bio_en):
    """Generate prerendered HTML for the bio section with both languages."""
    lines = []
    for text in bio_zh:
        lines.append(f'<p class="bio-paragraph zh">{esc(text)}</p>')
    for text in bio_en:
        lines.append(f'<p class="bio-paragraph en">{esc(text)}</p>')
    return '\n'.join(lines)


# ─── Image alt-text mapping for SEO ───

IMAGE_ALT_TEXT = {
    'images/honors/cute-alumni-group.webp': '2024 中國科技大學 59 週年校慶校友大會團體合照 — 黃崇源教授榮獲傑出校友',
    'images/honors/cute-distinguished-alumni.webp': '2024 中國科技大學傑出校友 — 黃崇源教授與校方師長合影',
    'images/honors/xtff-award.webp': '2024 臺北市杏壇芬芳獎頒獎典禮 — 黃崇源榮譽會長受獎',
    'images/honors/xtff-group.webp': '2024 臺北市杏壇芬芳獎頒獎典禮 — 黃崇源榮譽會長與家人師長合影',
    'images/honors/longshan-donation-1.webp': '2026 臺北市立龍山國中 — 黃崇源捐資母校受頒感謝狀現場照',
    'images/honors/longshan-donation-2.webp': '2026 臺北市立龍山國中 — 黃崇源捐資母校 115 學年感謝狀',
    'images/honors/laosong-130-1.webp': '2026 臺北市老松國民小學創校 130 周年校慶 — 黃崇源教授獲頒傑出校友貢獻獎榮譽狀（松小字第 1150005 號）',
    'images/honors/laosong-130-2.webp': '2026 臺北市老松國民小學創校 130 周年校慶 — 黃崇源教授榮獲傑出校友貢獻獎，與校長黎季昊合影',
    'images/honors/laosong-130-3.webp': '2026 臺北市老松國民小學創校 130 周年校慶 — 黃崇源先生獲頒金質獎捐款感謝狀（松小字第 1150006 號）',
    'images/honors/laosong-130-4.webp': '2026 臺北市老松國民小學創校 130 周年校慶 — 黃崇源先生榮獲金質獎捐款感謝狀，與校長黎季昊合影',
    'images/honors/laosong-advisor.webp': '2026 臺北市老松國民小學聘書 — 黃崇源先生榮膺 115 年校務顧問',
    'images/honors/laosong-honorary-president.webp': '2024 老松國小第 34 屆榮譽會長黃崇源教授參加校慶運動會',
    'images/honors/laosong-xlh-delegates.webp': '2023 臺北市國小學生家長會聯合會第 22 屆第一次代表大會合照',
    'images/honors/laosong-xlh-board.webp': '2023 臺北市國小學生家長會聯合會第一次理監事會議',
    'images/honors/laosong-xlh-election.webp': '2023 臺北市國小學生家長會聯合會第二十二屆常務理事選舉結果',
    'images/honors/laosong-xlh-west-district.webp': '2023 臺北市國小學生家長會聯合會第 22 屆西區理事選舉結果',
    'images/honors/laosong-president-award.webp': '2025 老松國小第 80 屆畢業典禮 — 黃崇源榮譽會長頒發榮譽會長獎',
    'images/honors/laosong-donation-1.webp': '2025 教師節 — 黃崇源榮譽會長捐贈新台幣伍拾萬元整予老松國小',
    'images/honors/laosong-donation-2.webp': '2025 教師節 — 老松國小校長與黃崇源榮譽會長合影留念',
    'images/honors/laosong-donation-3.webp': '2025 老松國小 114 學年度教師節慶祝活動海報',
    'images/honors/laosong-donation-4.webp': '2025 老松國小 114 年 9 月份捐款收支明細表',
    'images/honors/laosong-donation-5.webp': '2025 黃崇源榮譽會長捐贈老松國小新台幣伍拾萬元整支票',
    'images/honors/laosong-donation-6.webp': '2025 老松國小校長黎季昊頒贈黃崇源榮譽會長惠捐伍拾萬元感謝狀',
    'images/honors/solomon-newyear.webp': '2023 所羅門集團新年晚會 — 黃崇源教授擔任獨立董事',
    'images/honors/solomon-shareholders.webp': '2023 所羅門股份有限公司 112 年股東常會 — 黃崇源教授擔任獨立董事',
    'images/honors/ptphs-family.webp': '2006 中華民國斐陶斐榮譽學會獎 — 黃崇源教授與家人合影於國立交通大學頒獎典禮',
    'images/honors/ptphs-ceremony.webp': '2006 中華民國斐陶斐榮譽學會 95 年度榮譽會員授證典禮 — 黃崇源教授',
}


def extract_image_urls(html_str):
    """Extract all image URLs from showImageModal('...') calls in HTML.
    Strips ?query strings so alt-text lookups hit the dict key (e.g. `?v=2` cache-busters)."""
    urls = re.findall(r"showImageModal\('([^']+)'\)", html_str)
    return [u.split('?', 1)[0] for u in urls]


# ─── Parse object arrays from shared.js / stats-data.js ───

def parse_object_array(js_text, var_name, fields, source_label):
    """Parse a const array of flat objects, extracting the given string fields.
    Booleans named in `fields` with a trailing '?' are matched as `name: true`."""
    m = re.search(rf'const\s+{var_name}\s*=\s*\[', js_text)
    if not m:
        build_error(f"Could not find {var_name} in {source_label}")
        return []
    arr_start = m.end() - 1
    arr_end = find_bracket_end(js_text, arr_start)
    arr_content = js_text[arr_start + 1:arr_end - 1]

    items = []
    i = 0
    while i < len(arr_content):
        idx = arr_content.find('{', i)
        if idx == -1:
            break
        obj_end = find_bracket_end(arr_content, idx, '{', '}')
        obj_text = arr_content[idx:obj_end]
        i = obj_end

        item = {}
        for field in fields:
            if field.endswith('?'):
                name = field[:-1]
                item[name] = bool(re.search(rf'{name}:\s*true', obj_text))
            else:
                tm = re.search(rf'''{field}:\s*(['"])''', obj_text)
                item[field] = extract_js_string(obj_text, tm.end() - 1, tm.group(1))[0] if tm else None
        items.append(item)
    return items


# ─── Generate Stats HTML (#lab 統計卡預渲染) ───

def generate_stats_html(stats_zh, stats_en):
    """Replicate renderStats() markup for both languages so the 28 stat cards
    exist in static HTML for non-JS crawlers; JS re-renders them at runtime."""
    lines = []
    for lang, stats in (('zh', stats_zh), ('en', stats_en)):
        for it in stats:
            aria = f'View {it["label"]}' if lang == 'en' else f'點選查看 {it["label"]}'
            lines.append(
                f'<a class="stat-card {lang}" href="{esc(it["url"] or "#")}" target="_blank" '
                f'rel="noopener noreferrer" aria-label="{esc(aria)}">'
                f'<span class="stat-number">{esc(it["number"])}</span>'
                f'<div class="stat-label">{esc(it["label"])}</div></a>')
    return '\n'.join(lines)


# ─── Generate Nav HTML (無 JS 導覽後備) ───

def generate_nav_html(nav_zh, nav_en):
    """Static <a> navigation mirroring renderNav()'s skeleton (minus the JS-only
    mobile toggle). With JS on, renderNav() overwrites this at load; with JS off,
    the <noscript> style reveals all sections so the anchors work as jumps."""
    if len(nav_zh) != len(nav_en) or not nav_zh:
        build_error(f"navItemsZh ({len(nav_zh)}) / navItemsEn ({len(nav_en)}) mismatch or empty")
        return ''
    lines = ['<div class="container">',
             '    <div class="nav-container" id="nav-menu-container">',
             '        <ul class="nav-list">']
    for zh, en in zip(nav_zh, nav_en):
        if zh.get('section'):
            href, extra = '#' + zh['section'], ''
        else:
            href, extra = zh['href'], ' target="_blank" rel="noopener noreferrer"'
        label = f'<span class="zh">{esc(zh["label"])}</span><span class="en">{esc(en["label"])}</span>'
        lines.append(f'            <li><a class="nav-item" href="{esc(href)}"{extra}>{label}</a></li>')
    lines += ['        </ul>', '    </div>', '</div>']
    return '\n'.join(lines)


# ─── Generate Honors HTML ───

def generate_honors_html(honors_zh, honors_en):
    """Generate prerendered HTML for the honors list with both languages.
    Honors items contain raw HTML (buttons, spans) so we preserve them as-is.
    （原本額外塞入的隱藏 1x1 <img> 已移除：renderHonors() 執行時會將其刪除，
    會渲染 JS 的爬蟲反而看不到，且對螢幕報讀器是噪音；照片索引改由
    ImageGallery JSON-LD 與 sitemap 圖片條目完整涵蓋。）"""
    lines = []

    for html in honors_zh:
        lines.append(f'<li class="honors-item zh">{html}</li>')

    for html in honors_en:
        lines.append(f'<li class="honors-item en">{html}</li>')

    return '\n'.join(lines)


# ─── Generate Article HTML (Alumni / Honors recommendations) ───

def generate_article_html(paragraphs_zh, paragraphs_en):
    """Generate prerendered HTML for recommendation articles with both languages.
    Uses em-space indent (U+2003) matching the JS renderAlumniArticle/renderHonorsArticle."""
    lines = []
    for text in paragraphs_zh:
        lines.append(f'<p class="honors-paragraph zh">\u2003\u2003{esc(text)}</p>')
    for text in paragraphs_en:
        lines.append(f'<p class="honors-paragraph en">\u2003\u2003{esc(text)}</p>')
    return '\n'.join(lines)


# ─── Freshness stamps (--bump-dates) ───
# 同步 CLAUDE.md 慣例要求的五處日期戳記。採旗標而非每次建置自動更新，
# 以保留「同輸入必同輸出」的跨日冪等性；內容有實質變動時執行：
#   python3 build-prerender.py --bump-dates

def bump_dates(index_html):
    import datetime
    today = datetime.date.today().isoformat()

    index_html = re.sub(r'(<meta name="DCTERMS\.modified" content=")[0-9-]+(")',
                        rf'\g<1>{today}\g<2>', index_html)
    print(f"✅ DCTERMS.modified → {today}")

    file_patterns = [
        ('sitemap.xml',
         r'(<loc>https://canslab1\.github\.io/</loc>\s*<lastmod>)[0-9-]+(</lastmod>)',
         rf'\g<1>{today}\g<2>'),
        ('humans.txt', r'(Last update: )[0-9-]+', rf'\g<1>{today}'),
        ('llms-full.txt', r'(> Last updated: )[0-9-]+', rf'\g<1>{today}'),
        ('README.md', r'(\*\*最後更新\*\*：)[0-9-]+', rf'\g<1>{today}'),
    ]
    for fname, pattern, repl in file_patterns:
        path = os.path.join(BASE, fname)
        with open(path, 'r', encoding='utf-8') as f:
            content = f.read()
        new_content, n = re.subn(pattern, repl, content, count=1)
        if n:
            with open(path, 'w', encoding='utf-8') as f:
                f.write(new_content)
            print(f"✅ {fname} stamp → {today}")
        else:
            build_error(f"--bump-dates: stamp pattern not found in {fname}")

    return index_html


# ─── Main ───

def main():
    with open(os.path.join(BASE, 'papers-data.js'), 'r', encoding='utf-8') as f:
        papers_js = f.read()
    with open(os.path.join(BASE, 'projects-data.js'), 'r', encoding='utf-8') as f:
        projects_js = f.read()
    with open(os.path.join(BASE, 'shared.js'), 'r', encoding='utf-8') as f:
        shared_js = f.read()

    # Parse papers and projects
    papers_data = parse_papers_data(papers_js)
    projects_data = parse_projects_data(projects_js)

    total_papers = sum(len(item) for cat in papers_data for period in cat['periods'] for item in [period['items']])
    total_projects = sum(len(role['items']) for role in projects_data)
    print(f"Parsed: {total_papers} papers in {len(papers_data)} categories, {total_projects} projects in {len(projects_data)} roles")

    for i, cat in enumerate(papers_data):
        n = sum(len(p['items']) for p in cat['periods'])
        print(f"  Category {i+1}: {cat['titleZh']} — {n} papers, {len(cat['periods'])} periods")

    # Parse shared.js arrays
    print("\nParsing shared.js arrays...")
    bio_zh = parse_string_array(shared_js, 'bioZh')
    bio_en = parse_string_array(shared_js, 'bioEn')
    honors_zh = parse_string_array(shared_js, 'honorsZh')
    honors_en = parse_string_array(shared_js, 'honorsEn')
    alumni_zh = parse_string_array(shared_js, 'alumniArticleZh')
    alumni_en = parse_string_array(shared_js, 'alumniArticleEn')
    honors_article_zh = parse_string_array(shared_js, 'honorsArticleZh')
    honors_article_en = parse_string_array(shared_js, 'honorsArticleEn')
    nav_zh = parse_object_array(shared_js, 'navItemsZh', ['label', 'href', 'section', 'embed?'], 'shared.js')
    nav_en = parse_object_array(shared_js, 'navItemsEn', ['label', 'href', 'section', 'embed?'], 'shared.js')

    with open(os.path.join(BASE, 'stats-data.js'), 'r', encoding='utf-8') as f:
        stats_js = f.read()
    stats_zh = parse_object_array(stats_js, 'statsZh', ['number', 'label', 'url'], 'stats-data.js')
    stats_en = parse_object_array(stats_js, 'statsEn', ['number', 'label', 'url'], 'stats-data.js')
    print(f"  navItems: {len(nav_zh)}+{len(nav_en)}, stats: {len(stats_zh)}+{len(stats_en)} cards")
    print(f"  bioZh: {len(bio_zh)} paragraphs, bioEn: {len(bio_en)} paragraphs")
    print(f"  honorsZh: {len(honors_zh)} items, honorsEn: {len(honors_en)} items")
    print(f"  alumniArticleZh: {len(alumni_zh)} paragraphs, alumniArticleEn: {len(alumni_en)} paragraphs")
    print(f"  honorsArticleZh: {len(honors_article_zh)} paragraphs, honorsArticleEn: {len(honors_article_en)} paragraphs")

    # Generate HTML
    papers_html = generate_papers_html(papers_data)
    projects_html = generate_projects_html(projects_data)
    bio_html = generate_bio_html(bio_zh, bio_en)
    honors_html = generate_honors_html(honors_zh, honors_en)
    alumni_html = generate_article_html(alumni_zh, alumni_en)
    honors_article_html = generate_article_html(honors_article_zh, honors_article_en)
    stats_html = generate_stats_html(stats_zh, stats_en)
    nav_html = generate_nav_html(nav_zh, nav_en)

    jsonld = generate_jsonld(papers_data, projects_data)
    jsonld_str = json.dumps(jsonld, ensure_ascii=False, indent=4)

    index_path = os.path.join(BASE, 'index.html')
    with open(index_path, 'r', encoding='utf-8') as f:
        index_html = f.read()

    # Replace containers (idempotent — works whether each container is empty
    # or already filled by a previous run)
    containers = [
        ('<div id="papers-container">', papers_html, 'div', 'Papers'),
        ('<div id="projects-container">', projects_html, 'div', 'Projects'),
        ('<div id="bio-content">', bio_html, 'div', 'Bio'),
        ('<ul id="honors-list" class="honors-list">', honors_html, 'ul', 'Honors list'),
        ('<div id="alumni-article">', alumni_html, 'div', 'Alumni article'),
        ('<div id="honors-article">', honors_article_html, 'div', 'Honors article'),
        ('<div class="stats-grid" id="stats-container">', stats_html, 'div', 'Stats'),
        ('<nav id="main-nav">', nav_html, 'nav', 'Nav'),
    ]
    print()
    for open_tag, inner_html, tag_name, label in containers:
        index_html, ok = replace_container_inner(
            index_html, open_tag, '\n' + inner_html + '\n            ', tag_name)
        if ok:
            print('✅ %s container replaced' % label)
        else:
            build_error('%s container NOT found in index.html' % label)

    # Add JSON-LD before </body> — remove old auto-generated block first to prevent duplicates
    marker = '"Publications and Research Projects of Prof. Chung-Yuan Huang"'
    # Remove any existing auto-generated JSON-LD blocks
    while marker in index_html:
        start = index_html.rfind('<script type="application/ld+json">', 0, index_html.rfind(marker))
        end = index_html.find('</script>', index_html.rfind(marker)) + len('</script>')
        if start >= 0 and end > start:
            index_html = index_html[:start].rstrip() + '\n' + index_html[end:].lstrip()
        else:
            break
    ld_script = f'\n    <script type="application/ld+json">\n    {jsonld_str}\n    </script>'
    index_html = index_html.replace('</body>', ld_script + '\n</body>')
    print("✅ JSON-LD structured data added")

    n_articles = len([i for i in jsonld['itemListElement'] if i['item']['@type'] == 'ScholarlyArticle'])
    n_rp = len([i for i in jsonld['itemListElement'] if i['item']['@type'] == 'ResearchProject'])

    # 門檻斷言：解析結果遠低於已知規模即視為解析故障，而非真實內容變動
    if total_papers < 100:
        build_error(f"Sanity check failed: only {total_papers} papers parsed (expected 100+)")
    if total_projects < 30:
        build_error(f"Sanity check failed: only {total_projects} projects parsed (expected 30+)")
    if n_articles < 40:
        build_error(f"Sanity check failed: only {n_articles} ScholarlyArticle entries (expected 40+)")
    if not bio_zh or not bio_en or not honors_zh or not honors_en:
        build_error("Sanity check failed: bio/honors arrays are empty")
    if len(stats_zh) < 20 or len(stats_zh) != len(stats_en):
        build_error(f"Sanity check failed: stats cards {len(stats_zh)}/{len(stats_en)} (expected 20+, equal)")

    # 寫檔前的最終防線：任何錯誤都不得覆蓋 index.html
    if BUILD_ERRORS:
        print(f"\n❌ Build FAILED with {len(BUILD_ERRORS)} error(s); index.html NOT written:")
        for e in BUILD_ERRORS:
            print(f"  - {e}")
        sys.exit(1)

    if '--bump-dates' in sys.argv:
        index_html = bump_dates(index_html)
        if BUILD_ERRORS:
            print("\n❌ --bump-dates failed; index.html NOT written")
            sys.exit(1)

    with open(index_path, 'w', encoding='utf-8') as f:
        f.write(index_html)
    print(f"\n=== Summary ===")
    print(f"Pre-rendered papers: {total_papers}")
    print(f"Pre-rendered projects: {total_projects}")
    print(f"Pre-rendered bio: {len(bio_zh)}+{len(bio_en)} paragraphs")
    print(f"Pre-rendered honors: {len(honors_zh)}+{len(honors_en)} items")
    print(f"Pre-rendered alumni article: {len(alumni_zh)}+{len(alumni_en)} paragraphs")
    print(f"Pre-rendered honors article: {len(honors_article_zh)}+{len(honors_article_en)} paragraphs")
    print(f"JSON-LD: {n_articles} scholarly articles + {n_rp} research projects")
    print(f"index.html updated successfully!")

    print()
    minify_assets()


# ─── Minify deployed JS/CSS ───
# Sources stay readable (and build-prerender.py keeps parsing them);
# the HTML pages reference the generated *.min.* copies.

MINIFY_JS = ['shared.js', 'stats-data.js', 'papers-data.js', 'projects-data.js', 'stories.js']
MINIFY_CSS = ['shared.css', 'stories.css']


def minify_assets():
    try:
        import rjsmin
        import rcssmin
    except ImportError:
        rjsmin = rcssmin = None
        print("⚠️  rjsmin/rcssmin not installed (pip3 install --user rjsmin rcssmin);"
              " copying sources unminified")
    for name in MINIFY_JS + MINIFY_CSS:
        src = os.path.join(BASE, name)
        stem, ext = os.path.splitext(name)
        dst = os.path.join(BASE, stem + '.min' + ext)
        with open(src, 'r', encoding='utf-8') as f:
            code = f.read()
        if ext == '.js':
            out = rjsmin.jsmin(code) if rjsmin else code
        else:
            out = rcssmin.cssmin(code) if rcssmin else code
        with open(dst, 'w', encoding='utf-8') as f:
            f.write(out)
        print(f"✅ Minified {name}: {len(code):,} → {len(out):,} chars")


if __name__ == '__main__':
    main()
