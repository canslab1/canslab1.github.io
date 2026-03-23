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

用法：python3 build-prerender.py
"""

import json
import os
import re
import sys

BASE = os.path.dirname(os.path.abspath(__file__))


# ─────────────────────────────────────────────
# Parse JavaScript data using a bracket-based approach
# ─────────────────────────────────────────────

def extract_js_string(text, start):
    """Extract a JS single-quoted string starting at position `start` (must be the opening quote).
    Returns (string_value, end_position_after_closing_quote)."""
    assert text[start] == "'", f"Expected ' at pos {start}, got {text[start]!r}"
    i = start + 1
    chars = []
    while i < len(text):
        c = text[i]
        if c == '\\' and i + 1 < len(text):
            chars.append(text[i + 1])
            i += 2
        elif c == "'":
            return ''.join(chars), i + 1
        else:
            chars.append(c)
            i += 1
    raise ValueError("Unterminated string")


def find_bracket_end(text, start, open_b='[', close_b=']'):
    """Find the matching closing bracket, skipping strings."""
    assert text[start] == open_b
    depth = 1
    i = start + 1
    while i < len(text) and depth > 0:
        c = text[i]
        if c == "'":
            _, i = extract_js_string(text, i)
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
        print("ERROR: Could not find papersData")
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
        print("ERROR: Could not find projectsData")
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

    html = re.sub(r'\(([^()]*(?:SCI|SSCI|SCIE|EI)[^()]*)\)\s*\.?\s*$', r'<span class="paper-index">(\1)</span>', html)

    return html


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

    if papers_data:
        journal_cat = papers_data[0]
        for period in journal_cat['periods']:
            for text in period['items']:
                year_match = re.search(r'\((\d{4})\)', text)
                year = year_match.group(1) if year_match else None

                title_match = re.search(
                    r'\(\d{4}\)\s+(.+?)\.\s+(?:IEEE|PLoS|Physica|Scientific|Computational|Journal|Simulation|'
                    r'International|Applied|BioMed|Artificial|Discrete|WSEAS)',
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

    research_projects = []
    if projects_data:
        pi_cat = projects_data[0]
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
        print(f"  WARNING: Could not find {var_name}")
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


# ─── Generate Honors HTML ───

def generate_honors_html(honors_zh, honors_en):
    """Generate prerendered HTML for the honors list with both languages.
    Honors items contain raw HTML (buttons, spans) so we preserve them as-is."""
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

    jsonld = generate_jsonld(papers_data, projects_data)
    jsonld_str = json.dumps(jsonld, ensure_ascii=False, indent=4)

    index_path = os.path.join(BASE, 'index.html')
    with open(index_path, 'r', encoding='utf-8') as f:
        index_html = f.read()

    # Replace papers container
    old_papers = '<div id="papers-container"></div>'
    new_papers = '<div id="papers-container">\n' + papers_html + '\n            </div>'
    if old_papers in index_html:
        index_html = index_html.replace(old_papers, new_papers)
        print("\n✅ Papers container replaced")
    else:
        print("\n⚠️  Papers container not found (already replaced?)")

    # Replace projects container
    old_projects = '<div id="projects-container"></div>'
    new_projects = '<div id="projects-container">\n' + projects_html + '\n            </div>'
    if old_projects in index_html:
        index_html = index_html.replace(old_projects, new_projects)
        print("✅ Projects container replaced")
    else:
        print("⚠️  Projects container not found (already replaced?)")

    # Replace bio container
    old_bio = '<div id="bio-content"></div>'
    new_bio = '<div id="bio-content">\n' + bio_html + '\n            </div>'
    if old_bio in index_html:
        index_html = index_html.replace(old_bio, new_bio)
        print("✅ Bio container replaced")
    else:
        print("⚠️  Bio container not found (already replaced?)")

    # Replace honors list container
    old_honors = '<ul id="honors-list" class="honors-list"></ul>'
    new_honors = '<ul id="honors-list" class="honors-list">\n' + honors_html + '\n            </ul>'
    if old_honors in index_html:
        index_html = index_html.replace(old_honors, new_honors)
        print("✅ Honors list container replaced")
    else:
        print("⚠️  Honors list container not found (already replaced?)")

    # Replace alumni article container
    old_alumni = '<div id="alumni-article"></div>'
    new_alumni = '<div id="alumni-article">\n' + alumni_html + '\n            </div>'
    if old_alumni in index_html:
        index_html = index_html.replace(old_alumni, new_alumni)
        print("✅ Alumni article container replaced")
    else:
        print("⚠️  Alumni article container not found (already replaced?)")

    # Replace honors article container
    old_honors_art = '<div id="honors-article"></div>'
    new_honors_art = '<div id="honors-article">\n' + honors_article_html + '\n            </div>'
    if old_honors_art in index_html:
        index_html = index_html.replace(old_honors_art, new_honors_art)
        print("✅ Honors article container replaced")
    else:
        print("⚠️  Honors article container not found (already replaced?)")

    # Add JSON-LD before </body>
    ld_script = f'\n    <script type="application/ld+json">\n    {jsonld_str}\n    </script>'
    index_html = index_html.replace('</body>', ld_script + '\n</body>')
    print("✅ JSON-LD structured data added")

    with open(index_path, 'w', encoding='utf-8') as f:
        f.write(index_html)

    n_articles = len([i for i in jsonld['itemListElement'] if i['item']['@type'] == 'ScholarlyArticle'])
    n_rp = len([i for i in jsonld['itemListElement'] if i['item']['@type'] == 'ResearchProject'])
    print(f"\n=== Summary ===")
    print(f"Pre-rendered papers: {total_papers}")
    print(f"Pre-rendered projects: {total_projects}")
    print(f"Pre-rendered bio: {len(bio_zh)}+{len(bio_en)} paragraphs")
    print(f"Pre-rendered honors: {len(honors_zh)}+{len(honors_en)} items")
    print(f"Pre-rendered alumni article: {len(alumni_zh)}+{len(alumni_en)} paragraphs")
    print(f"Pre-rendered honors article: {len(honors_article_zh)}+{len(honors_article_en)} paragraphs")
    print(f"JSON-LD: {n_articles} scholarly articles + {n_rp} research projects")
    print(f"index.html updated successfully!")


if __name__ == '__main__':
    main()
