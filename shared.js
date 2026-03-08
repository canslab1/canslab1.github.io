/* ===== Academic Stats Data ===== */

const statsZh = [
    { number: '43篇',   label: '期刊論文', url: 'https://pure.lib.cgu.edu.tw/zh/persons/chung-yuan-huang-2/publications/?type=%2Fdk%2Fatira%2Fpure%2Fresearchoutput%2Fresearchoutputtypes%2Fcontributiontojournal%2Farticle' },
    { number: '49篇',   label: '國際研討會', url: 'https://pure.lib.cgu.edu.tw/zh/persons/chung-yuan-huang-2/publications/?type=%2Fdk%2Fatira%2Fpure%2Fresearchoutput%2Fresearchoutputtypes%2Fcontributiontoconference%2Finternational_conference_report' },
    { number: '10篇',   label: '專書專章', url: 'https://canslab1.github.io/CV.pdf' },
    { number: '19篇',   label: '國內研討會', url: 'https://canslab1.github.io/CV.pdf' },
    { number: '49次',   label: '文章及採訪', url: 'https://canslab1.github.io/CV.pdf' },
    { number: '22次',   label: '計畫主持人', url: 'https://pure.lib.cgu.edu.tw/zh/persons/chung-yuan-huang-2/projects/' },
    { number: '1408萬', label: '國科會預算', url: 'https://arspb.nstc.gov.tw/NSCWebFront/modules/talentSearch/talentSearch.do?action=initRsm17new&rsNo=2a2ef24adef742f58349c3c533bb7402&LANG=chi' },
    { number: '14次',   label: '共同主持人', url: 'https://canslab1.github.io/CV.pdf' },
    { number: '970次',  label: 'Citations', url: 'https://scholar.google.com/citations?user=0klfzfAAAAAJ&hl=en' },
    { number: '18篇',   label: 'h-index', url: 'https://scholar.google.com/citations?user=0klfzfAAAAAJ&hl=en' },
    { number: '30篇',   label: 'i10-index', url: 'https://scholar.google.com/citations?user=0klfzfAAAAAJ&hl=en' },
    { number: '51人',   label: '指導專題生', url: 'https://scholar.google.com/citations?user=0klfzfAAAAAJ&hl=en' },
    { number: '25人',   label: '指導碩士', url: 'https://scholar.google.com/citations?user=0klfzfAAAAAJ&hl=en' },
    { number: '2人',    label: '指導博士', url: 'https://scholar.google.com/citations?user=0klfzfAAAAAJ&hl=en' },
    { number: '68篇',   label: '論文審稿', url: 'https://scholar.google.com/citations?user=0klfzfAAAAAJ&hl=en' },
    { number: '27次',   label: '會議審稿', url: 'https://scholar.google.com/citations?user=0klfzfAAAAAJ&hl=en' },
    { number: '77次',   label: '議程委員', url: 'https://scholar.google.com/citations?user=0klfzfAAAAAJ&hl=en' },
    { number: '11次',   label: '計劃審查', url: 'https://scholar.google.com/citations?user=0klfzfAAAAAJ&hl=en' },
    { number: '10次',   label: '學術審查', url: 'https://scholar.google.com/citations?user=0klfzfAAAAAJ&hl=en' },
    { number: '35次',   label: '專題演講', url: 'https://scholar.google.com/citations?user=0klfzfAAAAAJ&hl=en' },
    { number: '10年',   label: '教授', url: 'https://scholar.google.com/citations?user=0klfzfAAAAAJ&hl=en' },
    { number: '5年',    label: '副教授', url: 'https://scholar.google.com/citations?user=0klfzfAAAAAJ&hl=en' },
    { number: '5年',    label: '助理教授', url: 'https://scholar.google.com/citations?user=0klfzfAAAAAJ&hl=en' },
    { number: '5年',    label: '講師', url: 'https://scholar.google.com/citations?user=0klfzfAAAAAJ&hl=en' },
    { number: '11年',   label: '業界經驗', url: 'https://scholar.google.com/citations?user=0klfzfAAAAAJ&hl=en' },
    { number: '4個',    label: '合聘單位', url: 'https://canslab1.github.io/CV.pdf' },
    { number: '3年',    label: '獨立董事', url: 'https://canslab1.github.io/CV.pdf' },
    { number: '7年',    label: '學生家長會', url: 'https://canslab1.github.io/CV.pdf' }
];

const statsEn = [
    { number: '43',   label: 'Journal Papers', url: 'https://canslab1.github.io/CV.pdf#page=21' },
    { number: '49',   label: "Int'l Conferences", url: 'https://canslab1.github.io/CV.pdf#page=24' },
    { number: '10',   label: 'Book Chapters', url: 'https://canslab1.github.io/CV.pdf#page=27' },
    { number: '19',   label: 'Dom. Conferences', url: 'https://canslab1.github.io/CV.pdf#page=29' },
    { number: '49',   label: 'Commentaries & Op-Eds', url: 'https://scholar.google.com/citations?user=0klfzfAAAAAJ&hl=en' },
    { number: '22',   label: 'PI Projects', url: 'https://canslab1.github.io/CV.pdf#page=3' },
    { number: '$14M', label: 'NSTC Budget', url: 'https://canslab1.github.io/CV.pdf#page=3' },
    { number: '14',   label: 'Co-PI Projects', url: 'https://canslab1.github.io/CV.pdf#page=4' },
    { number: '970',  label: 'Citations', url: 'https://scholar.google.com/citations?user=0klfzfAAAAAJ&hl=en' },
    { number: '18',   label: 'h-index', url: 'https://scholar.google.com/citations?user=0klfzfAAAAAJ&hl=en' },
    { number: '30',   label: 'i10-index', url: 'https://scholar.google.com/citations?user=0klfzfAAAAAJ&hl=en' },
    { number: '51',   label: 'UG Students', url: 'https://scholar.google.com/citations?user=0klfzfAAAAAJ&hl=en' },
    { number: '25',   label: 'Msc Students', url: 'https://scholar.google.com/citations?user=0klfzfAAAAAJ&hl=en' },
    { number: '2',    label: 'PhD Students', url: 'https://scholar.google.com/citations?user=0klfzfAAAAAJ&hl=en' },
    { number: '68',   label: 'Paper Reviews', url: 'https://scholar.google.com/citations?user=0klfzfAAAAAJ&hl=en' },
    { number: '27',   label: 'Conference Reviews', url: 'https://scholar.google.com/citations?user=0klfzfAAAAAJ&hl=en' },
    { number: '77',   label: 'Program Committee', url: 'https://scholar.google.com/citations?user=0klfzfAAAAAJ&hl=en' },
    { number: '11',   label: 'Project Reviews', url: 'https://scholar.google.com/citations?user=0klfzfAAAAAJ&hl=en' },
    { number: '10',   label: 'Academic Reviews', url: 'https://scholar.google.com/citations?user=0klfzfAAAAAJ&hl=en' },
    { number: '35',   label: 'Invited Talks', url: 'https://scholar.google.com/citations?user=0klfzfAAAAAJ&hl=en' },
    { number: '10 yrs',   label: 'Full Professor', url: 'https://scholar.google.com/citations?user=0klfzfAAAAAJ&hl=en' },
    { number: '5 yrs',    label: 'Assoc. Prof.', url: 'https://scholar.google.com/citations?user=0klfzfAAAAAJ&hl=en' },
    { number: '5 yrs',    label: 'Asst. Prof.', url: 'https://scholar.google.com/citations?user=0klfzfAAAAAJ&hl=en' },
    { number: '5 yrs',    label: 'Lecturer', url: 'https://scholar.google.com/citations?user=0klfzfAAAAAJ&hl=en' },
    { number: '11 yrs',   label: 'Industry Experience', url: 'https://scholar.google.com/citations?user=0klfzfAAAAAJ&hl=en' },
    { number: '4',    label: 'Joint Appointments', url: 'https://canslab1.github.io/CV.pdf#page=1' },
    { number: '3 yrs',    label: 'Independent Director', url: 'https://canslab1.github.io/CV.pdf#page=1' },
    { number: '7 yrs',    label: 'Parent Association', url: 'https://canslab1.github.io/CV.pdf#page=1' }
];

/* ===== Honors Data ===== */

const honorsZh = [
    '<span style="color: #2980b9;">中國科技大學</span><br>2024.11 榮獲傑出校友<br>2000.11 榮獲優秀校友<br>校史唯一先後獲頒傑出與優秀校友',
    '<span style="color: #2980b9;">臺北市杏壇芬芳獎</span><br>2024.11 榮獲臺北市教育局公開表揚',
    '<span style="color: #2980b9;">臺北市老松國民小學學生家長會</span><br>2024.10 榮膺第 34 屆榮譽會長<br>2023.10 擔任臺北市小聯會常務理事<br>2022.10 擔任第 32 和 33 屆會長<br>113 學年起提供畢業生榮譽會長獎，鼓勵多元學習表現優異學童',
    '<span style="color: #2980b9;">所羅門股份有限公司</span><br>2022.06 榮膺第十二屆董事會獨立董事<br>2022.06 擔任薪酬、審計、永續發展委員',
    '<span style="color: #2980b9;">中華民國斐陶斐榮譽學會獎</span><br>2006.06 榮獲國立交通大學公開頒發<br>表彰卓越學術成就與優異品德操守',
    '<span style="color: #2980b9;">臺北市杏壇芬芳獎推薦文</span><br>源源不絕無私大愛，守護美好優質老松'
];

const honorsEn = [
    '<span style="color: #2980b9;">China University of Technology</span><br>2024.11 Distinguished Alumni Award<br>2000.11 Outstanding Alumni Award<br>the only recipient in university history of both Distinguished and Outstanding Alumni honors.',
    '<span style="color: #2980b9;">2024.10 Taipei Honored Contributors to Education</span><br>Bottomless Devotion to Taipei Lao-Song Elementary School',
    '<span style="color: #2980b9;">Laosong Elementary School, Taipei City</span><br>2024.10 Honorary President, 34th Parents&#39; Association<br>2022.10-2024.10 President of the 32nd and 33rd Parents&#39; Association<br>Starting from the 113th Academic Year, the Honorary PTA President\'s Award will be presented to graduating students',
    '<span style="color: #2980b9;">Solomon Technology Corporation</span><br>2022.06 Served as Independent Director on the 12th Board of Directors of Solomon Co., Ltd.<br>2022.06 Served on the Compensation, Audit, and Sustainability Committees of the 12th Board of Directors of Solomon Co., Ltd.',
    '<span style="color: #2980b9;">Phi Tau Phi Honor Award</span><br>2006.06 Awarded by National Chiao Tung University,<br>in recognition of academic excellence and exemplary character',
    '<span style="color: #2980b9;">Recommendation Article for the Taipei Honored contributors to education</span>'
];

/* ===== Navigation Data ===== */

const navItemsZh = [
    { label: '介紹', href: null, section: 'overview' },
    { label: '實驗室', href: 'https://canslab1.github.io/lab.html', page: 'lab' },
    { label: '著作', href: 'https://scholar.google.com/citations?user=0klfzfAAAAAJ&hl=en' },
    { label: '計畫', href: 'https://arspb.nstc.gov.tw/NSCWebFront/modules/talentSearch/talentSearch.do?action=initRsm17new&rsNo=2a2ef24adef742f58349c3c533bb7402&LANG=chi' },
    { label: '履歷', href: 'https://canslab1.github.io/CV.pdf' },
    { label: '維基', href: 'https://sites.google.com/view/gscott-huang' },
    { label: '臉書', href: 'https://www.facebook.com/gscott.huang/' },
    { label: '故事', href: 'https://canslab1.github.io/stories.html' },
    { label: '相簿', href: 'https://www.dropbox.com/scl/fo/ejwke56gkn1erv7meid4k/AHCPQXAw_RpiEn_PeiBfCuI?rlkey=b8k879fufdh98x2s3fzkatk5y&e=1&dl=0' },
    { label: '評論', href: 'https://www.google.com/search?q=%22Huang+Chung-yuan%22+%22Taipei+Times%22+%22Chang+Gung%22+site%3Awww.taipeitimes.com' },
    { label: '投書', href: 'https://www.google.com/search?q=%22%E9%BB%83%E5%B4%87%E6%BA%90%22+%22%E8%87%AA%E7%94%B1%E6%99%82%E5%A0%B1%22+site%3Atalk.ltn.com.tw' }
];

const navItemsEn = [
    { label: 'About', href: null, section: 'overview' },
    { label: 'Lab', href: 'https://canslab1.github.io/lab.html', page: 'lab' },
    { label: 'Papers', href: 'https://scholar.google.com/citations?user=0klfzfAAAAAJ&hl=en' },
    { label: 'Projects', href: 'https://arspb.nstc.gov.tw/NSCWebFront/modules/talentSearch/talentSearch.do?action=initRsm17new&rsNo=2a2ef24adef742f58349c3c533bb7402&LANG=chi' },
    { label: 'CV', href: 'https://canslab1.github.io/CV.pdf' },
    { label: 'Wiki', href: 'https://sites.google.com/view/gscott-huang' },
    { label: 'FB', href: 'https://www.facebook.com/gscott.huang/' },
    { label: 'Stories', href: 'https://canslab1.github.io/stories.html' },
    { label: 'Photos', href: 'https://www.dropbox.com/scl/fo/ejwke56gkn1erv7meid4k/AHCPQXAw_RpiEn_PeiBfCuI?rlkey=b8k879fufdh98x2s3fzkatk5y&e=1&dl=0' },
    { label: 'Press', href: 'https://www.google.com/search?q=%22Huang+Chung-yuan%22+%22Taipei+Times%22+%22Chang+Gung%22+site%3Awww.taipeitimes.com' },
    { label: 'Op-Eds', href: 'https://www.google.com/search?q=%22%E9%BB%83%E5%B4%87%E6%BA%90%22+%22%E8%87%AA%E7%94%B1%E6%99%82%E5%A0%B1%22+site%3Atalk.ltn.com.tw' }
];

/* ===== Footer Data ===== */

const footerLinks = [
    { href: 'https://www.cgu.edu.tw/', label: '前往長庚大學官方網站', src: 'https://canslab1.github.io/images/cgu.png', alt: '長庚大學 logo' },
    { href: 'https://www.cgu.edu.tw/csie/', label: '前往長庚大學資訊工程學系', src: 'https://canslab1.github.io/images/csie.png', alt: '長庚資工系 logo' },
    { href: 'https://web.tlsps.tp.edu.tw/nss/p/index', label: '前往老松國小', src: 'https://canslab1.github.io/images/laosong.png', alt: '老松國小 logo' },
    { href: 'https://www.cute.edu.tw/', label: '前往中國科技大學', src: 'https://canslab1.github.io/images/cute.jpg', alt: '中國科技大學 logo' },
    { href: 'https://xingtan.tiec.tp.edu.tw/Register/Profile/20240715152009039506', label: '前往臺北市杏壇芬芳錄', src: 'https://canslab1.github.io/images/favicon.ico', alt: '臺北市杏壇芬芳錄 logo' },
    { href: 'https://orcid.org/0000-0002-8680-6755', label: '前往ORCiD', src: 'https://orcid.org/assets/vectors/orcid.logo.icon.svg', alt: 'ORCiD logo' }
];

/* ===== Language Toggle ===== */

function toggleLanguage() {
    const html = document.documentElement;
    const btn = document.querySelector('.language-toggle');
    const isZh = html.lang === 'zh-TW';

    html.lang = isZh ? 'en' : 'zh-TW';
    btn.textContent = isZh ? '中文' : 'English';

    document.querySelectorAll('.en').forEach(el => {
        el.style.display = isZh ? '' : 'none';
    });
    document.querySelectorAll('.zh').forEach(el => {
        el.style.display = isZh ? 'none' : '';
    });

    renderStats();
    renderHonors();
    renderNav();
}

/* ===== Section Toggle ===== */

function showSection(sectionId, event) {
    document.querySelectorAll('.section').forEach(s => s.classList.remove('active'));
    document.querySelectorAll('.nav-item').forEach(n => n.classList.remove('active'));

    const target = document.getElementById(sectionId);
    if (target) target.classList.add('active', 'fade-in');
    if (event && event.target) event.target.classList.add('active');

    const navContainer = document.querySelector('.nav-container');
    if (navContainer && navContainer.classList.contains('show')) {
        navContainer.classList.remove('show');
    }

    window.scrollTo({ top: 0, behavior: 'smooth' });
}

function toggleNav() {
    document.querySelector('.nav-container').classList.toggle('show');
}

/* ===== Render Functions ===== */

function renderStats() {
    const isEn = document.documentElement.lang === 'en';
    const stats = isEn ? statsEn : statsZh;
    const container = document.getElementById('stats-container');
    if (!container) return;
    container.innerHTML = '';

    stats.forEach(item => {
        const card = document.createElement('button');
        card.className = 'stat-card';
        card.setAttribute('aria-label', isEn ? `View ${item.label}` : `點選查看 ${item.label}`);
        card.innerHTML = `<span class="stat-number">${item.number}</span><div class="stat-label">${item.label}</div>`;
        card.addEventListener('click', () => window.open(item.url, '_blank'));
        container.appendChild(card);
    });
}

function renderHonors() {
    const isEn = document.documentElement.lang === 'en';
    const honors = isEn ? honorsEn : honorsZh;
    const container = document.getElementById('honors-list');
    if (!container) return;
    container.innerHTML = '';

    honors.forEach(txt => {
        const li = document.createElement('li');
        li.innerHTML = txt;
        container.appendChild(li);
    });
}

/* currentPage: 'index' or 'lab' */
let _currentPage = 'index';

function renderNav() {
    const nav = document.getElementById('main-nav');
    if (!nav) return;

    const isEn = document.documentElement.lang === 'en';
    const items = isEn ? navItemsEn : navItemsZh;

    nav.innerHTML = `
        <div class="container">
            <button class="nav-toggle" onclick="toggleNav()">
                <span class="zh">☰ 導覽</span><span class="en">☰ Navigation</span>
            </button>
            <div class="nav-container"></div>
        </div>
    `;

    // Fix language display inside nav-toggle
    const toggleSpans = nav.querySelectorAll('.nav-toggle span');
    toggleSpans.forEach(span => {
        if (span.classList.contains('en')) {
            span.style.display = isEn ? 'inline' : 'none';
        } else {
            span.style.display = isEn ? 'none' : 'inline';
        }
    });

    const container = nav.querySelector('.nav-container');

    items.forEach((item, i) => {
        if (item.section) {
            // "介紹/Overview" button — links to index page section
            if (_currentPage === 'index') {
                const btn = document.createElement('button');
                btn.className = 'nav-item active';
                btn.textContent = item.label;
                btn.addEventListener('click', (e) => showSection(item.section, e));
                container.appendChild(btn);
            } else {
                const a = document.createElement('a');
                a.className = 'nav-item';
                a.href = 'https://canslab1.github.io/';
                a.textContent = item.label;
                container.appendChild(a);
            }
        } else if (item.page === 'lab') {
            if (_currentPage === 'lab') {
                const a = document.createElement('a');
                a.className = 'nav-item active';
                a.href = '#';
                a.setAttribute('aria-current', 'page');
                a.textContent = item.label;
                container.appendChild(a);
            } else {
                const a = document.createElement('a');
                a.className = 'nav-item';
                a.href = item.href;
                a.target = '_blank';
                a.rel = 'noopener noreferrer';
                a.textContent = item.label;
                container.appendChild(a);
            }
        } else {
            const a = document.createElement('a');
            a.className = 'nav-item';
            a.href = item.href;
            a.target = '_blank';
            a.rel = 'noopener noreferrer';
            a.textContent = item.label;
            container.appendChild(a);
        }
    });
}

function renderFooter() {
    const footer = document.getElementById('main-footer');
    if (!footer) return;

    const div = document.createElement('div');
    div.className = 'container footer-content';

    footerLinks.forEach(link => {
        const a = document.createElement('a');
        a.href = link.href;
        a.target = '_blank';
        a.rel = 'noopener noreferrer';
        a.setAttribute('aria-label', link.label);

        const img = document.createElement('img');
        img.src = link.src;
        img.alt = link.alt;
        img.className = 'footer-logo';

        a.appendChild(img);
        div.appendChild(a);
    });

    footer.appendChild(div);
}

/* ===== Initialization ===== */

function initShared(currentPage) {
    _currentPage = currentPage || 'index';
    renderStats();
    renderHonors();
    renderNav();
    renderFooter();
}

/* ===== Microsoft Clarity ===== */

(function(c,l,a,r,i,t,y){
    c[a]=c[a]||function(){(c[a].q=c[a].q||[]).push(arguments)};
    t=l.createElement(r);t.async=1;t.src="https://www.clarity.ms/tag/"+i+"?ref=bwt";
    y=l.getElementsByTagName(r)[0];y.parentNode.insertBefore(t,y);
})(window, document, "clarity", "script", "rzlnthqbys");
