/* ===== Honors Data ===== */

const honorsZh = [
    '<span class="honors-highlight">中國科技大學</span><br>2024.11 榮獲傑出校友 <button class="interview-btn" onclick="showVideoModal(\'028zexOcXZo\')">專訪</button> <button class="interview-btn" onclick="showImageModal(\'images/honors/cute-alumni-group.webp\')">照片一</button> <button class="interview-btn" onclick="showImageModal(\'images/honors/cute-distinguished-alumni.webp\')">照片二</button><br>2000.11 榮獲優秀校友<br>校史唯一先後獲頒傑出與優秀校友',
    '<span class="honors-highlight">臺北市杏壇芬芳獎</span><br>2024.11 榮獲臺北市教育局公開表揚 <button class="interview-btn" onclick="showImageModal(\'images/honors/xtff-award.webp\')">照片一</button> <button class="interview-btn" onclick="showImageModal(\'images/honors/xtff-group.webp\')">照片二</button><br>源源不絕無私大愛，守護美好優質老松',
    '<span class="honors-highlight">臺北市老松國民小學</span><br>2026.03 榮膺 115 年校務顧問 <button class="interview-btn" onclick="showImageModal(\'images/honors/laosong-advisor.webp\')">照片一</button><br>2025.10 教師節捐款新台幣伍拾萬元整 <button class="interview-btn" onclick="showImageModal(\'images/honors/laosong-donation-1.webp\')">照片一</button> <button class="interview-btn" onclick="showImageModal(\'images/honors/laosong-donation-2.webp\')">照片二</button> <button class="interview-btn" onclick="showImageModal(\'images/honors/laosong-donation-3.webp\')">照片三</button> <button class="interview-btn" onclick="showImageModal(\'images/honors/laosong-donation-4.webp\')">照片四</button> <button class="interview-btn" onclick="showImageModal(\'images/honors/laosong-donation-5.webp\')">照片五</button><br>111 學年起，連續四學年捐款總額超過新臺幣壹佰萬元整',
    '<span class="honors-highlight">臺北市老松國民小學學生家長會</span><br>2024.10 榮膺第 34 屆榮譽會長 <button class="interview-btn" onclick="showImageModal(\'images/honors/laosong-honorary-president.webp\')">照片</button><br>2023.10 擔任臺北市小聯會常務理事 <button class="interview-btn" onclick="showImageModal(\'images/honors/laosong-xlh-delegates.webp\')">照片一</button> <button class="interview-btn" onclick="showImageModal(\'images/honors/laosong-xlh-board.webp\')">照片二</button> <button class="interview-btn" onclick="showImageModal(\'images/honors/laosong-xlh-election.webp\')">照片三</button> <button class="interview-btn" onclick="showImageModal(\'images/honors/laosong-xlh-west-district.webp\')">照片四</button><br>2022.10 擔任第 32 和 33 屆會長<br>113 學年起設立畢業生榮譽會長獎，每班壹名，鼓勵多元學習表現優異學童 <button class="interview-btn" onclick="showImageModal(\'images/honors/laosong-president-award.webp\')">照片</button>',
    '<span class="honors-highlight">所羅門股份有限公司</span><br>2022.06 榮膺第十二屆董事會獨立董事 <button class="interview-btn" onclick="showImageModal(\'images/honors/solomon-newyear.webp\')">照片一</button> <button class="interview-btn" onclick="showImageModal(\'images/honors/solomon-shareholders.webp\')">照片二</button><br>2022.06 擔任薪酬、審計、永續發展委員',
    '<span class="honors-highlight">中華民國斐陶斐榮譽學會獎</span><br>2006.06 榮獲國立交通大學公開頒發 <button class="interview-btn" onclick="showImageModal(\'images/honors/ptphs-family.webp\')">照片一</button> <button class="interview-btn" onclick="showImageModal(\'images/honors/ptphs-ceremony.webp\')">照片二</button><br>表彰卓越學術成就與優異品德操守'
];

const honorsEn = [
    '<span class="honors-highlight">China University of Technology</span><br>2024.11 Distinguished Alumni Award <button class="interview-btn" onclick="showVideoModal(\'028zexOcXZo\')">Interview</button> <button class="interview-btn" onclick="showImageModal(\'images/honors/cute-alumni-group.webp\')">Photo 1</button> <button class="interview-btn" onclick="showImageModal(\'images/honors/cute-distinguished-alumni.webp\')">Photo 2</button><br>2000.11 Outstanding Alumni Award<br>the only recipient in university history of both Distinguished and Outstanding Alumni honors.',
    '<span class="honors-highlight">2024.10 Taipei Honored Contributors to Education</span><br>2024.11 Publicly commended by the Taipei City Department of Education <button class="interview-btn" onclick="showImageModal(\'images/honors/xtff-award.webp\')">Photo 1</button> <button class="interview-btn" onclick="showImageModal(\'images/honors/xtff-group.webp\')">Photo 2</button><br>Bottomless Devotion to Taipei Lao-Song Elementary School',
    '<span class="honors-highlight">Laosong Elementary School, Taipei City</span><br>2026.03 Appointed as School Affairs Advisor for the year 2026 <button class="interview-btn" onclick="showImageModal(\'images/honors/laosong-advisor.webp\')">Photo 1</button><br>2025.10 Donated NT$500,000 on Teachers&#39; Day <button class="interview-btn" onclick="showImageModal(\'images/honors/laosong-donation-1.webp\')">Photo 1</button> <button class="interview-btn" onclick="showImageModal(\'images/honors/laosong-donation-2.webp\')">Photo 2</button> <button class="interview-btn" onclick="showImageModal(\'images/honors/laosong-donation-3.webp\')">Photo 3</button> <button class="interview-btn" onclick="showImageModal(\'images/honors/laosong-donation-4.webp\')">Photo 4</button> <button class="interview-btn" onclick="showImageModal(\'images/honors/laosong-donation-5.webp\')">Photo 5</button><br>Since the 111th Academic Year, cumulative donations over four consecutive years exceeded NT$1,000,000',
    '<span class="honors-highlight">Laosong Elementary School Parents&#39; Association, Taipei City</span><br>2024.10 Honorary President, 34th Parents&#39; Association <button class="interview-btn" onclick="showImageModal(\'images/honors/laosong-honorary-president.webp\')">Photo</button><br>2023.10 Standing Director of Taipei Elementary School PTA Federation <button class="interview-btn" onclick="showImageModal(\'images/honors/laosong-xlh-delegates.webp\')">Photo 1</button> <button class="interview-btn" onclick="showImageModal(\'images/honors/laosong-xlh-board.webp\')">Photo 2</button> <button class="interview-btn" onclick="showImageModal(\'images/honors/laosong-xlh-election.webp\')">Photo 3</button> <button class="interview-btn" onclick="showImageModal(\'images/honors/laosong-xlh-west-district.webp\')">Photo 4</button><br>2022.10-2024.10 President of the 32nd and 33rd Parents&#39; Association<br>Starting from the 113th Academic Year, the Honorary PTA President\'s Award will be presented to graduating students <button class="interview-btn" onclick="showImageModal(\'images/honors/laosong-president-award.webp\')">Photo</button>',
    '<span class="honors-highlight">Solomon Technology Corporation</span><br>2022.06 Served as Independent Director on the 12th Board of Directors of Solomon Co., Ltd. <button class="interview-btn" onclick="showImageModal(\'images/honors/solomon-newyear.webp\')">Photo 1</button> <button class="interview-btn" onclick="showImageModal(\'images/honors/solomon-shareholders.webp\')">Photo 2</button><br>2022.06 Served on the Compensation, Audit, and Sustainability Committees of the 12th Board of Directors of Solomon Co., Ltd.',
    '<span class="honors-highlight">Phi Tau Phi Honor Award</span><br>2006.06 Awarded by National Chiao Tung University, <button class="interview-btn" onclick="showImageModal(\'images/honors/ptphs-family.webp\')">Photo 1</button> <button class="interview-btn" onclick="showImageModal(\'images/honors/ptphs-ceremony.webp\')">Photo 2</button><br>in recognition of academic excellence and exemplary character'
];

/* ===== Navigation Data ===== */

const navItemsZh = [
    { label: '介紹', href: null, section: 'bio' },
    { label: '榮譽', href: null, section: 'overview' },
    { label: '綜覽', href: null, section: 'lab' },
    { label: '著作', href: null, section: 'papers' },
    { label: '計畫', href: null, section: 'projects' },
    { label: '程式', href: null, section: 'software' },
    { label: '文章', href: null, section: 'articles' },
    { label: '評論', href: null, section: 'press' },
    { label: '履歷', href: 'https://canslab1.github.io/CV.pdf', embed: true },
    { label: '故事', href: 'https://canslab1.github.io/stories.html' },
    { label: '相簿', href: 'https://www.dropbox.com/scl/fo/ejwke56gkn1erv7meid4k/AHCPQXAw_RpiEn_PeiBfCuI?rlkey=b8k879fufdh98x2s3fzkatk5y&e=1&dl=0' }
];

const navItemsEn = [
    { label: 'About', href: null, section: 'bio' },
    { label: 'Honors', href: null, section: 'overview' },
    { label: 'Overview', href: null, section: 'lab' },
    { label: 'Papers', href: null, section: 'papers' },
    { label: 'Projects', href: null, section: 'projects' },
    { label: 'Software', href: null, section: 'software' },
    { label: 'Articles', href: null, section: 'articles' },
    { label: 'Op-Eds', href: null, section: 'press' },
    { label: 'CV', href: 'https://canslab1.github.io/CV.pdf', embed: true },
    { label: 'Stories', href: 'https://canslab1.github.io/stories.html' },
    { label: 'Photos', href: 'https://www.dropbox.com/scl/fo/ejwke56gkn1erv7meid4k/AHCPQXAw_RpiEn_PeiBfCuI?rlkey=b8k879fufdh98x2s3fzkatk5y&e=1&dl=0' }
];

/* ===== Footer Data ===== */

const footerLinks = [
    { href: 'https://www.cgu.edu.tw/csie/', label: '前往長庚大學資訊工程學系', labelEn: 'CGU CSIE', src: 'https://canslab1.github.io/images/csie-thumb.png', alt: '長庚資工系 logo' },
    { href: 'https://web.tlsps.tp.edu.tw/nss/p/index', label: '前往老松國小', labelEn: 'Laosong Elementary', src: 'https://canslab1.github.io/images/laosong-thumb.png', alt: '老松國小 logo' },
    { href: 'https://xingtan.tiec.tp.edu.tw/Register/Profile/20240715152009039506', label: '前往臺北市杏壇芬芳錄', labelEn: 'Taipei Honored Contributors', src: 'https://canslab1.github.io/images/favicon-thumb.ico', alt: '臺北市杏壇芬芳錄 logo' },
    { href: 'https://orcid.org/0000-0002-8680-6755', label: '前往ORCiD', labelEn: 'ORCiD', src: 'https://orcid.org/assets/vectors/orcid.logo.icon.svg', alt: 'ORCiD logo' },
    { href: 'https://www.facebook.com/gscott.huang/', label: '前往 Facebook', labelEn: 'Facebook', src: 'https://canslab1.github.io/images/facebook-thumb.png', alt: 'Facebook logo' },
    { href: 'https://sites.google.com/view/gscott-huang', label: '前往 Google Sites', labelEn: 'Google Sites', src: 'https://canslab1.github.io/images/google-thumb.png', alt: 'Google logo' },
    { href: 'https://scholar.google.com/citations?user=0klfzfAAAAAJ&hl=en', label: '前往 Google Scholar', labelEn: 'Google Scholar', src: 'https://canslab1.github.io/images/google-scholar-thumb.png', alt: 'Google Scholar logo' },
    { href: 'https://github.com/canslab1', label: '前往 GitHub', labelEn: 'GitHub', src: 'https://canslab1.github.io/images/github-thumb.png', alt: 'GitHub logo' },
    { href: 'https://pure.lib.cgu.edu.tw/zh/persons/chung-yuan-huang-2', label: '前往長庚大學 Pure 學術檔案', labelEn: 'CGU Pure', src: 'https://canslab1.github.io/images/cgu-thumb.png', alt: '長庚大學 Pure logo' }
];

/* ===== Language Toggle ===== */

function toggleLanguage() {
    const html = document.documentElement;
    const btn = document.querySelector('.language-toggle');
    const isZh = html.lang === 'zh-TW';

    html.lang = isZh ? 'en' : 'zh-TW';
    btn.textContent = isZh ? '中文' : 'English';

    renderStats();
    renderHonors();
    renderAlumniArticle();
    renderHonorsArticle();
    renderBio();
    renderPapers();
    renderProjects();
    renderNav();
    renderFooter();
    updateAriaLabels();
}

/* ===== Section Toggle ===== */

function showSection(sectionId, event) {
    document.querySelectorAll('.section').forEach(s => s.classList.remove('active', 'fade-in'));
    document.querySelectorAll('.nav-item').forEach(n => n.classList.remove('active'));

    const target = document.getElementById(sectionId);
    if (target) target.classList.add('active', 'fade-in');
    if (event && event.target) event.target.classList.add('active');

    _activeSection = sectionId;

    /* 切回原有區段時清空 iframe 以釋放資源 */
    const frame = document.getElementById('embed-frame');
    if (frame) frame.src = 'about:blank';

    /* 重設文章 PDF 檢視器 */
    hideArticlePdf();

    const navContainer = document.querySelector('.nav-container');
    if (navContainer && navContainer.classList.contains('show')) {
        navContainer.classList.remove('show');
    }

    window.scrollTo({ top: 0, behavior: 'smooth' });
}

function showEmbed(url, event) {
    document.querySelectorAll('.section').forEach(s => s.classList.remove('active', 'fade-in'));
    document.querySelectorAll('.nav-item').forEach(n => n.classList.remove('active'));

    const embed = document.getElementById('embed');
    if (embed) embed.classList.add('active', 'fade-in');

    const frame = document.getElementById('embed-frame');
    if (frame) frame.src = url;

    if (event && event.target) event.target.classList.add('active');

    _activeSection = 'embed';

    const navContainer = document.querySelector('.nav-container');
    if (navContainer && navContainer.classList.contains('show')) {
        navContainer.classList.remove('show');
    }

    window.scrollTo({ top: 0, behavior: 'smooth' });
}

function toggleNav() {
    const navContainer = document.querySelector('.nav-container');
    if (navContainer) navContainer.classList.toggle('show');
}

/* ===== ARIA Labels ===== */

function updateAriaLabels() {
    const isEn = document.documentElement.lang === 'en';

    const lab = document.getElementById('lab');
    if (lab) lab.setAttribute('aria-label', isEn ? 'Overview' : '綜覽');

    const embed = document.getElementById('embed');
    if (embed) embed.setAttribute('aria-label', isEn ? 'Embedded content' : '嵌入內容');
}

/* ===== Render Functions ===== */

function renderStats() {
    const container = document.getElementById('stats-container');
    if (!container) return;
    const isEn = document.documentElement.lang === 'en';
    const stats = isEn ? statsEn : statsZh;
    container.innerHTML = '';

    stats.forEach(item => {
        const card = document.createElement('a');
        card.className = 'stat-card';
        card.href = item.url || '#';
        card.target = '_blank';
        card.rel = 'noopener noreferrer';
        card.setAttribute('aria-label', isEn ? `View ${item.label}` : `點選查看 ${item.label}`);
        card.innerHTML = `<span class="stat-number">${item.number}</span><div class="stat-label">${item.label}</div>`;
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

var _activeSection = 'bio';
function renderNav() {
    const nav = document.getElementById('main-nav');
    if (!nav) return;

    const isEn = document.documentElement.lang === 'en';
    const items = isEn ? navItemsEn : navItemsZh;

    nav.innerHTML = `
        <div class="container">
            <button class="nav-toggle">
                <span class="zh">☰ 導覽</span><span class="en">☰ Navigation</span>
            </button>
            <div class="nav-container"></div>
        </div>
    `;

    const navToggleBtn = nav.querySelector('.nav-toggle');
    if (navToggleBtn) navToggleBtn.addEventListener('click', toggleNav);

    const container = nav.querySelector('.nav-container');

    const ul = document.createElement('ul');
    ul.className = 'nav-list';
    ul.setAttribute('role', 'menubar');

    items.forEach(item => {
        const li = document.createElement('li');
        li.setAttribute('role', 'none');

        if (item.section) {
            const btn = document.createElement('button');
            btn.className = 'nav-item' + (item.section === _activeSection ? ' active' : '');
            btn.setAttribute('role', 'menuitem');
            btn.textContent = item.label;
            btn.addEventListener('click', (e) => showSection(item.section, e));
            li.appendChild(btn);
        } else if (item.embed) {
            const btn = document.createElement('button');
            btn.className = 'nav-item' + (_activeSection === 'embed' ? ' active' : '');
            btn.setAttribute('role', 'menuitem');
            btn.textContent = item.label;
            btn.addEventListener('click', (e) => showEmbed(item.href, e));
            li.appendChild(btn);
        } else {
            const a = document.createElement('a');
            a.className = 'nav-item';
            a.setAttribute('role', 'menuitem');
            a.href = item.href;
            a.target = 'canslab-' + item.label;
            a.rel = 'noopener noreferrer';
            a.textContent = item.label;
            li.appendChild(a);
        }

        ul.appendChild(li);
    });

    container.appendChild(ul);
}

function renderFooter() {
    const footer = document.getElementById('main-footer');
    if (!footer) return;

    footer.innerHTML = '';
    const isEn = document.documentElement.lang === 'en';

    const nav = document.createElement('nav');
    nav.setAttribute('aria-label', isEn ? 'Partner links' : '合作單位連結');

    const div = document.createElement('div');
    div.className = 'container footer-content';

    footerLinks.forEach(link => {
        const a = document.createElement('a');
        a.href = link.href;
        a.target = '_blank';
        a.rel = 'noopener noreferrer';
        a.className = 'footer-link';
        const ariaText = isEn ? link.labelEn : link.label;
        a.setAttribute('aria-label', ariaText);

        const tooltip = document.createElement('span');
        tooltip.className = 'footer-tooltip';
        tooltip.textContent = isEn ? link.labelEn : link.label.replace(/^前往\s?/, '');
        a.appendChild(tooltip);

        const img = document.createElement('img');
        img.src = link.src;
        img.alt = link.alt;
        img.className = 'footer-logo';

        a.appendChild(img);
        div.appendChild(a);
    });

    nav.appendChild(div);
    footer.appendChild(nav);
}

/* ===== Honors Article Data ===== */

const honorsArticleZh = [
    '「蘇○○，早安！何○○，早安！」',
    '如冬日暖陽般親切地喊出學童姓名的崇源會長，每天早晨總是站在穿堂旁學生家長會辦公室門口，以充滿朝氣的笑容，和藹可親的問候，迎接每一位踏入老松校園的學生。在學童的眼中，他是永遠面帶微笑，樂意彎下腰、蹲下身，傾聽他們說話，關懷他們心情的學生家長會長。',
    '崇源會長是老松校友，畢業於民國七十一年，家族三代都就讀老松，舅舅更長期擔任老松國小訓導主任。崇源會長常掛在嘴邊的一句話：「老松國小的大小事，都是我的大事，能夠協助母校親師生，是畢生的榮幸！」',
    '崇源會長帶著一顆感恩、奉獻、回饋的心，投入老松學生家長會志工服務的行列，期間經歷了家長會委員、志工總幹事、活動組長、常務委員、副會長等職務，非常熟悉會務。在疫情接近尾聲之際，崇源會長結合了一群來自低、中、高年級，橫跨園藝、保健、導護、愛閱、圖書、ＥＱ、學輔所有組別的熱忱家長志工，秉持著自發、互動、共好的教育理念，入校服務，全力支援老松國小於疫後陸續恢復的各項活動，積極協助學童在花木扶疏、綠意盎然的校園環境中探索與學習。',
    '時間是最寶貴的資源。為了協助疫後重返校園的老松師生與愛心志工，崇源會長毫無私心，奉獻出寶貴的時間與無比的愛心，每日坐鎮學生家長會辦公室，用充滿活力與正向積極的任事態度，激勵團隊夥伴，透過集思廣益，淬鍊出多元的趣味創意與靈活想法，盡心盡力辦好學生家長會務，協助校務與教育工作的發展，關心老松孩童的校園生活。',
    '曾是國小校警，如今擔任保全的蕭○○先生，常對崇源會長豎起大姆指說：「會長，您真的很厲害！整日守護在學校，處理會務，關懷親師生，其他學校的家長會辦公室通常大門深鎖很少使用，只有開會和活動才有人氣，兩者差異很大。」',
    '除了獻出寶貴的時間，崇源會長也以大學教授平日從事研究的嚴謹態度，檢視學生家長會辦理的親子成長、親職講座、志工成長、感恩餐會、一日遊等活動。從前置準備，到執行細節，再到預期結果，他總是不厭其煩和團隊夥伴反覆推敲活動的每項細節，可能的意外，以及臨場可行的補救之道。崇源會長常對團隊夥伴說：「舉辦活動，就要辦得熱鬧有意義，讓參與活動的親師生與愛心志工收獲滿滿，身心靈滿載而歸。」',
    '無論在公開場合，還是線上家長討論群組，家長對於高情商的崇源會長的真心讚美不知凡幾。大家的推崇包括以下但不止於：老松國小擁有崇源會長，是親師生的好福氣；崇源會長是一位回來老松國小造福的好會長；崇源會長親和力十足，沒有距離感，是值得安心信賴的會長；崇源會長是我們眼中最用心、最認真、親力親為的滿分會長；從來不曾有過像現在的感受，原來會長和家長的關係可以如此緊密，感謝崇源會長為老松親師生所做的一切，以及造成的美好改變。',
    '羅馬不是一天造成的。家長溢於言表的由衷讚美之辭，源自於崇源會長多年來表裡如一的個人特質。誠如他的團隊夥伴在長時間相處及第一手觀察之後，都異口同聲形容：「崇源會長為人謙沖有禮，出錢出力，不求回報。他是一位願意傾聽、樂於溝通、即時回應、詳細說明、熱心協助、將心比心、身先士卒、以身作則、無私奉獻的會長，是家長與學校之間通暢的溝通者。」陽光又虛懷若谷的崇源會長則是透過實際的作為，屢次在公開場合強調，他是老松國小家長的駐校代表，是家長派駐在老松國小服務親師生的頭號值日生，是可愛的老松學弟妹及家長安心可靠的好夥伴。',
    '在一些家長的眼中，崇源會長好比聞聲救苦，有大智慧的長者。疫後時期，雖然學習活動漸漸步上軌道，但是每位學童的恢復過程快慢不一，許多孩子的學習會突然卡關，家長面臨此類問題需要校方協助，同樣來自教育體系，在大學執教二十餘年的崇源會長總是站在學生學習的教育立場思考，不嫌麻煩主動和家長溝通，詳盡回覆，然後協助家長，向校方傳達問題緣由，提醒關鍵重點，委請學校處室負責人員適時伸出援手，給予學童和家長應有且適當的助益。',
    '不只是家長信賴的代言人與值日生，崇源會長也是學校政策和活動的最佳行銷者。為了促成家長和校方的教育理念、政策方向與教學步伐一致，崇源會長致力於擔任老松親師生之間最佳的溝通橋樑與潤滑劑，協助處室主任，透過適當的傳播管道，向家長詳盡地解釋校方近期正在推展的雙語教學等重要政策及可能的影響層面，以及即將舉辦的國際交流等重要活動及其背後的宗旨和教育目的。除此之外，崇源會長也主動定期向家長說明政策執行與活動辦理的精彩過程與豐碩成果，推展親師合力，共好共榮。',
    '美好的教育環境的營造，需要所有的教育夥伴同心協力，克盡職責，同時，教育現場需要典範，指引教育夥伴堅定向前。崇源會長一直是老松國小正向發展中最給力的教育合夥人，推薦崇源會長參加本年度杏壇芬芳錄評選實為老松國小的榮幸，亦是親師生的期待。',
    '老松國小學生家長會長黃崇源教授，是一位足以成為眾人楷模的教育合夥人，若能獲得委員的肯定與表彰，寫入杏壇芬芳錄裡，必能激發其他教育工作夥伴的熱忱與活力，禮讚杏壇，齊力實現美好健全的教育理想！'
];

const honorsArticleEn = [
    'Good morning, Su! Good morning, Ho!',
    'Every morning, President Chung-Yuan stands at the entrance of the Parents\u2019 Association office next to the auditorium, and greets every student entering the Lao-Song Campus by name and with a warm, cheerful smile. In the minds of Taipei Lao-Song Elementary School students, he is the man who always has a smile on his face, and who is always willing to bend down to the students\u2019 level to hear what\u2019s on their minds.',
    'President Chung-Yuan graduated from Lao-Song in 1982, one of three generations in his family to attend. His uncle has been the Director of Student Affairs at the school for many years. That\u2019s why President Chung-Yuan frequently says, \u201cEverything that happens at Lao-Song is my business, and it\u2019s my honor to help my alma mater, its teachers, and its students.\u201d',
    'His gratitude, dedication, and positive feedback is evident in his devotion to serving as a volunteer for the Lao-song Parents\u2019 Association. Over the years he has served as the Association\u2019s Vice-President, Standing Committee member, and Chief of Volunteer Activities and Volunteer Programs. He knows the business of the Association inside and out. Toward the end of the epidemic he came into contact with parents of children at all grade levels, with backgrounds in horticulture, health care, guidance counseling, academic counseling, and books and reading, and promoted the educational philosophy of self-initiative, interaction, and the common good. The goal: to actively assist children in exploring and learning.',
    'President Chung-Yuan has always been willing to share one of his most valuable resources: time. To help Lao-Song students, teachers, and volunteers with the transition back to in-class learning after the COVID-19 pandemic, he selflessly took on the task of sitting in the Parents\u2019 Association office every day, encouraging those around him with his energy, positive attitude, and skill at brainstorming creative ideas. Showing his care for the lives and efforts of Lao-Song students, he helped reorganize the Parents\u2019 Association to assist in all school affairs, especially educational work.',
    'Mr. Hsiao, a former police officer assigned to Lao-Song who now serves as a security guard, often points his finger at President Chung-Yuan and says, \u201cPresident, you are really great! You\u2019re at the school all day long, holding meetings and caring for teachers and students. In other schools the PA office is usually locked and only used for meetings and events. There\u2019s a big difference between the two!\u201d',
    'In addition to contributing his valuable time, President Chung-Yuan oversees parent-child development, parenting seminars, volunteer development, Thanksgiving meals, day trips, and other PA-sponsored activities with the same rigor that university professors dedicate to their research. From preparation to implementation, he takes great pain to review all activity details with his fellow volunteers, deal with surprises, and come up with remedies for problems. He often tells his teammates, \u201cWhen organizing an event, we must make it lively and meaningful, so that everyone who participates can feel the full benefits, and return home with a sense of fulfillment.\u201d',
    'Whether in person or in online discussion groups, parents have given President Chung-Yuan countless compliments for the way he deals with them and others. Some example comments include, \u201cLao-Song Elementary is blessed to have President Chung-Yuan\u201d; \u201cPresident Chung-Yuan is full of empathy, with no sense of separating differences, and can be trusted to achieve peace of mind\u201d; \u201che is the most sincere, serious, and personally involved president\u201d; \u201cthere has never been a relationship between the PA president and parents like the one we have right now\u201d; and \u201cwe thank him for the wonderful changes he has brought about.\u201d',
    'Rome was not built in a day, and the parents of Lao-Song students clearly appreciate President Chung-Yuan\u2019s long-term commitment to the school, in addition to his personal qualities. After spending time with him and observing him firsthand, some of his fellow volunteers described him as \u201chumble and courteous, contributing without asking for anything in return\u201d; \u201cwilling to listen, willing to communicate, responds immediately, explains in detail, is eager to help\u201d; \u201cputs his heart and soul into his work, leads by good example\u201d; and \u201ca skilled communicator who connects parents with teachers and school employees.\u201d He cheerfully describes himself in public as the Number One Duty Officer assigned by Lao-Song parents to serve teachers and students at the school, and a reliable partner for both parents and teachers.',
    'In the minds of some parents, President Chung-Yuan acts like a wise elder who listens and responds with the goal of reducing suffering. He understands that students in post-pandemic classrooms are recovering academically at different rates, with a significant number of children feeling stuck. He also understands that the parents of such students need help from the school, and maybe from other sources. As a university educator for over 20 years, he thinks in terms of the learning needs of students. He willfully takes the time to reply to the questions of parents in detail, as well as assist them in navigating the various offices of school authorities. He helps parents convey their concerns and key points to school administrators, and provides assistance to staff members if they request it.',
    'As a trusted parent advocate and duty officer, President Chung-Yuan is a skilled marketer of school policies and programs. To promote consistency between the needs of parents and the educational philosophy of Lao-Song, he is committed to serving as a safe bridge between teachers and students, assisting school officials with the task of explaining their policies to parents in detail, using appropriate communication channels that respect the needs of all parties. He is especially concerned about supporting the success of new initiatives such as bilingual education and international exchanges, explaining their objectives and possible impacts. In addition, he goes out of his way to keep parents informed of progress and the fruitful results of programs and policies, both old and new, thus promoting parent-teacher cooperation for the common good of all in the Lao-Song community.',
    'Creating positive educational environments requires cooperation among all concerned partners, as well as the presence of role models to guide forward progress. Outside of teachers and staff, President Chung-Yuan has been the most active and contributing agent to the positive development of Lao-Song Elementary. Therefore we, the teachers and students of Lao-Song, enthusiastically endorse his participation in this year\u2019s Taipei Honored Contributors to Education Selection. Both the education sector and Lao-Song Elementary will benefit from his involvement.'
];

/* ===== Alumni Article Data ===== */

const alumniArticleZh = [
    '余等謹此誠摯推薦民國71年第37屆畢業生黃崇源教授作為臺北市老松國小創校130週年校慶傑出校友貢獻獎之表揚人選。黃教授現任長庚大學資訊工程學系教授、所羅門股份有限公司獨立董事，長期致力於人工智慧領域之學術研究和人才培育，曾於2024年同時榮獲中國科技大學傑出校友及臺北市杏壇芬芳錄兩項肯定並公開表揚，近年返回母校擔任學生家長會會長、榮譽會長。時至今日，始終秉持食果子拜樹頭、食米飯拜田頭的初衷，感謝母校恩師的極力栽培和持續鼓勵，並竭力協助可愛的學弟妹及尊敬的師長，其事蹟可謂楷模，誠為本獎最具代表性的候選人。',
    '黃教授學術成就卓越，著作等身，歷年主持國科會與疾管署專題研究計畫逾三十件，主題涵蓋新興流行性感染症傳播動態建模、公共衛生政策評估與人工智慧前沿演算法，歷年成果廣受肯定，均發表於國際頂尖學術期刊，總引用次數破千。其研究亦曾榮獲模擬學會高被引論文獎及多項國際研討會最佳論文獎。尤其難能可貴的是，黃教授將先進電腦模擬路數應用於公衛議題，對台灣防疫政策與自動體外除顫器設置提供具體貢獻，彰顯其專業不限於學術象牙塔，更深植於社會實踐中。',
    '促使黃教授走上學術研究卓越生涯之啟蒙點，正是老松國小。在其發表於《老松兒童》的兩篇文章中，余等得知小學三年級導師蔡幸美老師不僅每日親授課後指導，更鼓勵他閱讀大量課外讀物，種下了快速閱讀與知識探索的種子；四年級的姚嘉奉老師鼓勵其創作現代詩，五、六年級的陳艦東與呂御榮老師則以生動活潑的教學法激發其對語文的熱情與無窮想像。這些經驗奠定其邏輯思維與語文運用能力，也讓他深刻感受到教育之力量，並堅信有為者亦若是。其今日於大學場域大放異彩的邏輯思辨、研究創造與跨域敘事等能力，無一不可追溯至老松校園裡那段恩師接力栽培的童年歲月。',
    '近年，黃教授除了持續在人工智慧領域深耕與貢獻，更全心全力投入老松國小學生家長會，歷任活動組長、志工總幹事、常務委員、副會長，乃至會長與榮譽會長，與數十位跨年級、跨組別的家長志工協同服務親師生，協助推動雙語教學、校園安全、親職講座與親子成長活動。他每日駐守會辦，親迎學生上下學、傾聽家長意見、協助教師推動校務發展，並於疫後積極復辦愛心志工旅遊、感恩餐會及敬師活動，充分展現高度參與及責任感。誠如臺北市杏壇芬芳錄推薦文所述，其將心比心、身先士卒精神屢獲校內外讚譽，家長與師長咸稱黃會長是老松國小最可靠的教育夥伴。',
    '努力兮！努力兮！學求優，行求良\u2014這段取自老松國小校歌的優美詞句，正是黃教授長年實踐的座右銘。他不僅在學術上追求卓越優秀，也在人才培育與社會服務中實踐「行求良」之崇高理想。黃教授始終認為，自己的每一分成就與社會肯定，都是小學時期恩師的啟發與鼓舞所致。如今，他用自己的力量一點一滴回饋母校，守護母校，也用具體的行動成為眾多老松學童日後仰望的典範。',
    '綜上所述，黃教授的學術貢獻、社會服務與母校情誼三者交織成一幅老松國小130週年校慶傑出校友貢獻獎的最佳圖像。誠盼委員明鑑，惠予肯定，使其精神風範得以為後學仿效，並為我老松130週年校慶添光增彩！'
];

const alumniArticleEn = [
    'We enthusiastically endorse Professor Chung-Yuan Huang (\u9ec3\u5d07\u6e90\u6559\u6388) as a candidate for a Distinguished Alumni Contribution Award on the occasion of Lao-Song\u2019s 130th anniversary.',
    'A graduate of our school\u2019s 37th cohort (1982), Professor Huang has received two similar awards. In 2024 he was honored as both Distinguished Alumnus (\u5091\u51fa\u6821\u53cb) of China University of Technology (\u4e2d\u570b\u79d1\u6280\u5927\u5b78), and a recipient of Taipei Honored Contributors to Education (\u674f\u58c7\u82ac\u82b3\u9304). He is currently a Professor in the Department of Computer Science and Information Engineering at Chang Gung University, and an Independent Director of the Solomon Technology Corporation. He has a long history of academic research, education, and social engagement in the field of artificial intelligence. Yet despite all of his responsibilities, he made time to serve as a president and honorary president of the Lao-Song Parent\u2019s Association. In that manner he continues to express his gratitude to the nurturing teachers and mentors who supported him, abiding by the principle, \u201cWhen eating fruit, thank the tree; when eating rice, thank the field.\u201d',
    'Professor Huang has served as principal investigator for more than thirty research projects funded by Taiwan\u2019s National Science and Technology Council and Centers for Disease Control. His studies address topics such as the transmission dynamics of emerging infectious diseases, public health policy evaluation, and frontier artificial intelligence algorithms. His results have consistently been published in top-tier academic journals, with more than 1,000 citations. He received a High Impact Paper Award from the Society for Modeling and Simulation, and multiple Best Paper Awards from international conferences. Of particular note is his application of advanced computational simulation methods to public health topics, thus contributing to government policies for epidemic prevention and automated external defibrillator deployment. These accomplishments reflect his motivation for social engagement in his scholarly work.',
    'The origins of his academic career took place at Lao-Song Elementary. In an article he wrote for Lao-Song Children, he describes the encouragement he received from his third grade teacher, Ms. Sing-Mei Tsai (\u8521\u5e78\u7f8e\u8001\u5e2b), who not only provided him with daily tutoring, but also encouraged him to read widely beyond the classroom, sowing the seeds of his curiosity and his ability to quickly absorb large amounts of written information. His fourth grade teacher, Mr. Jiia-Fong Yao (\u59da\u5609\u5949\u8001\u5e2b), encouraged him to write poetry, and the vivid and lively teaching methods of his fifth and sixth grade teachers, Mr. Jian-Dong Chen (\u9673\u8266\u6771\u8001\u5e2b) and Mr. Yuh-Rong Lyu (\u5442\u5fa1\u69ae\u8001\u5e2b), sparked his love for language. In short, the logical rigor, research creativity, and interdisciplinary fluency that define his scholarship can be traced to his dedicated teachers during his years as a Lao-Song student.',
    'During the past seven years Professor Huang far exceeded his usual volunteer time in various Lao-Song Parent Association positions. He successively took on duties as activity director, chief volunteer coordinator, executive committee member, vice president, president, and honorary president. He collaborated with dozens of parents to promote bilingual education, campus safety campaigns, parenting seminars, and family development programs. He personally greeted students at arrival, and said goodbye at dismissal. He became a daily fixture in the school office, listening to the concerns of parents, assisting teachers with administrative tasks, and renewing volunteer appreciation banquets and teacher recognition events after the pandemic. His Taipei Honored Contributors to Education commendation cited his empathy and willingness to lead by example.',
    'Professor Huang\u2019s support embodies the ideal of moral integrity through educational practice and civic engagement. He continues to give back to his alma mater through tangible actions that offer a strong and positive role model for Lao-Song students. As we observe the 130th anniversary of the school\u2019s distinctive educational work, we encourage the committee to recognize his contributions so that his example can continue to serve as an inspiration.'
];

const honorsFigures = [
    {
        src: 'https://canslab1.github.io/images/honors1.webp',
        alt: '老松國小．學生家長會．黃崇源會長過往二年的部分成果展現',
        captionZh: '▲ 在老松學童的眼中，崇源會長就是聖誕老公公，帶來歡笑，帶來令人難忘的回憶。',
        captionEn: '▲ To the children of Lao-Song Elementary, President Chung-Yuan is their Santa Claus\u2014bringing joy, spreading laughter, and leaving behind memories that linger long after the moment has passed.'
    },
    {
        src: 'https://canslab1.github.io/images/honors2.webp',
        alt: '老松國小．學生家長會．黃崇源會長過往二年的部分成果展現',
        captionZh: '▲ 崇源會長出錢出力，盡心盡力辦好學生家長會主辦的親師生及愛心志工活動。',
        captionEn: '▲ President Chung-Yuan contributed both time and resources with wholehearted dedication, ensuring that every event for teachers, students, and volunteers was carried out with care, creativity, and purpose.'
    },
    {
        src: 'https://canslab1.github.io/images/honors3.webp',
        alt: '老松國小．學生家長會．黃崇源會長過往二年的部分成果展現',
        captionZh: '▲ 在崇源會長付出寶貴的時間與無比的愛心後，老松國小親師生同心協力共好共榮。',
        captionEn: '▲ Through the generous time and boundless care of President Chung-Yuan, the Lao-Song community\u2014teachers, parents, and students alike\u2014came together in harmony, working hand in hand toward shared growth and collective flourishing.'
    },
    {
        src: 'https://canslab1.github.io/images/honors4.webp',
        alt: '老松國小．學生家長會．黃崇源會長過往二年的部分成果展現',
        captionZh: '▲ 崇源會長永遠面帶微笑，和老松所有學童（包含幼兒園幼童）相處融洽。',
        captionEn: '▲ With a warm and ever-present smile, President Chung-Yuan connects effortlessly with every child at Lao-Song\u2014from kindergarteners to upper-grade students\u2014creating bonds built on kindness, trust, and joy.'
    }
];

/* ===== Biography Data ===== */

const bioZh = [
    '黃崇源博士（1970年8月13日生）擔任多項職務，包括大學教授、資訊科學家、醫學中心研究員及上市公司獨立董事。黃博士目前居住在臺北市萬華區，擔任長庚大學資訊工程學系教授，同時兼任人工智慧學程與研究所、人工智慧研究中心合聘教授等職務。',
    '1985年，他考入中國工商專校電子資料處理科，在全職工作的同時，歷經七年半工半讀完成學士學位。畢業後，他在澎湖防衛司令部服義務兵役，然後在宇博電腦（臺灣）股份有限公司擔任四年的軟體工程師。1998年，他考取新竹國立交通大學計算機資訊科學所，並於2000年與2005年分別獲頒碩士和博士學位。在攻讀博士學位期間，擔任中國科技大學資訊工程學系專任講師。2006年，他被授予中華民國斐陶斐學會學術榮譽會員資格。獲得博士學位後，黃崇源在元培醫科大學資訊工程學系擔任一學年的專任助理教授，並兼任系主任和生物醫學資訊技術學群的副執行長。自2006年起，黃博士擔任長庚大學資訊工程學系專任助理教授，兼任資訊中心教學服務組組長。2010年升等副教授，2015年升等正教授。自2019年8月起，他在長庚紀念醫院神經內科擔任聯合聘任研究員。',
    '黃博士目前的研究興趣包括複雜系統計算建模與模擬、社會建模與模擬及流行病傳播動力學。他的研究結合了遺傳演算法、學習分類系統、複雜網絡分析等人工智慧方面的最新研究成果。他的第一篇學術論文（關於SARS傳播動態與公衛政策模擬）發表在2004年的學術期刊JASSS上。此後，他發表了170篇學術論文，其中44篇被SCI/SSCI收錄，49篇發表在國際會議論文集上，10篇發表在書籍章節中，67篇發表在國內會議論文和大眾報紙專欄文章中。根據Google Scholar的統計，截至2026年2月，他發表的論文被引用994次（其中30篇每篇被引用10次以上），h-index為18。'
];

const bioEn = [
    'The many roles played by Chung-Yuan Huang, Ph.D. (b. August 13, 1970) include university professor, information scientist, medical center research fellow, and independent director of a publicly listed company. Currently residing in the Wanhua District of Taipei, Dr. Huang is a full professor in the Department of Computer Science and Information Engineering at Chang Gung University, concurrently holding a joint appointment as professor in the Artificial Intelligence undergraduate program and research institute, and the Artificial Intelligence Research Center.',
    'In 1985 he was admitted to the Department of Electronic Data Processing at China Junior College of Technology and Commerce, completing his bachelor\u2019s degree over seven years of part-time study while working full-time. After graduating he completed two years of mandatory military service at the Penghu Defense Command Headquarters, and then worked for four years (1994-98) as a software engineer at IPACS Computer and Service (Taiwan) Ltd. In 1998 he was accepted into the Department of Computer Information Science of National Chiao Tung University in Hsinchu, where he earned his master\u2019s (2000) and doctoral degrees (2005). During his doctoral studies he served as a full-time lecturer in the Department of Computer Science and Information Engineering at China University of Technology. In 2006 he was awarded membership in the Republic of China Phi Tau Phi Scholastic Honor Society. After receiving his Ph.D., Chung-Yuan spent one year as a full-time assistant professor in the Department of Computer Science and Information Engineering at Yuanpei University of Medical Technology, also serving as department head and Deputy Chief Executive Officer of the school\u2019s Biomedical Information Technology Cluster. Since 2006, Dr. Huang has been a full-time assistant professor in the Department of Computer Science and Information Engineering at Chang Gung University, and head of the university\u2019s Computer Center Instructional Services Group. He was promoted to associate professor in 2010 and full professor in 2015. Since August 2019 he has fulfilled a joint appointment as a researcher in the Department of Neurology at Chang Gung Memorial Hospital.',
    'Dr. Huang\u2019s current research interests include computational modeling and simulation of complex systems, social modeling and simulation, and epidemic transmission dynamics. His research increasingly incorporates the latest findings in genetic algorithms, learning classifier systems, complex network analysis, and other aspects of artificial intelligence. His first academic paper (on the effectiveness of SARS containment policies) appeared in a 2004 edition of the Journal of Artificial Societies and Social Simulation. Since then, he has published 170 research papers, including 44 indexed by SCI/SSCI, 49 in international conference proceedings, 10 book chapters, and 67 domestic conference papers and mass-circulation newspaper op-eds. According to Google Scholar, as of February 2026 his publications have been cited 994 times (30 of them more than ten times each), with an h-index of 18.'
];

function renderBio() {
    const container = document.getElementById('bio-content');
    if (!container) return;
    const isEn = document.documentElement.lang === 'en';
    const paragraphs = isEn ? bioEn : bioZh;
    container.innerHTML = '';
    paragraphs.forEach(text => {
        const p = document.createElement('p');
        p.className = 'bio-paragraph';
        p.textContent = text;
        container.appendChild(p);
    });
}

function renderAlumniArticle() {
    const container = document.getElementById('alumni-article');
    if (!container) return;

    const isEn = document.documentElement.lang === 'en';
    const paragraphs = isEn ? alumniArticleEn : alumniArticleZh;

    container.innerHTML = '';

    paragraphs.forEach(text => {
        const p = document.createElement('p');
        p.className = 'honors-paragraph';
        p.textContent = '\u2003\u2003' + text;
        container.appendChild(p);
    });
}

function renderHonorsArticle() {
    const container = document.getElementById('honors-article');
    if (!container) return;

    const isEn = document.documentElement.lang === 'en';
    const paragraphs = isEn ? honorsArticleEn : honorsArticleZh;

    container.innerHTML = '';

    paragraphs.forEach(text => {
        const p = document.createElement('p');
        p.className = 'honors-paragraph';
        p.textContent = '\u2003\u2003' + text;
        container.appendChild(p);
    });

    honorsFigures.forEach(fig => {
        const figure = document.createElement('figure');

        const img = document.createElement('img');
        img.src = fig.src;
        img.alt = fig.alt;
        img.className = 'honors-image';
        img.loading = 'lazy';
        figure.appendChild(img);

        const caption = document.createElement('figcaption');
        caption.className = 'honors-caption';

        const zhSpan = document.createElement('span');
        zhSpan.className = 'zh';
        zhSpan.textContent = fig.captionZh;
        caption.appendChild(zhSpan);

        const enSpan = document.createElement('span');
        enSpan.className = 'en';
        enSpan.textContent = fig.captionEn;
        caption.appendChild(enSpan);

        figure.appendChild(caption);
        container.appendChild(figure);
    });
}

/* ===== Papers Rendering ===== */

function formatPaperText(text) {
    var html = text
        .replace(/Huang C\.Y\.(\*?)/g, '<strong>Huang C.Y.$1</strong>')
        .replace(/黃崇源/g, '<strong>黃崇源</strong>');

    html = html.replace(/(https?:\/\/[^\s,)]+[^\s,).])/g, '<a href="$1" target="_blank" rel="noopener noreferrer">$1</a>');
    html = html.replace(/doi:(10\.\d{4,}\/[^\s,)]+[^\s,).])/g, '<a href="https://doi.org/$1" target="_blank" rel="noopener noreferrer">doi:$1</a>');
    html = html.replace(/\[PDF:([^\]]+)\]/g, '<a href="$1" target="_blank" rel="noopener noreferrer">[PDF]</a>');

    html = html.replace(/\(([^()]*(?:SCI|SSCI|SCIE|EI)[^()]*)\)\s*\.?\s*$/, '<span class="paper-index">($1)</span>');

    return html;
}

function renderPapers() {
    var container = document.getElementById('papers-container');
    if (!container || typeof papersData === 'undefined') return;

    var isEn = document.documentElement.lang === 'en';
    container.innerHTML = '';

    papersData.forEach(function(category) {
        var catDiv = document.createElement('div');
        catDiv.className = 'papers-category';

        var catTitle = document.createElement('h3');
        catTitle.className = 'papers-category-title';
        catTitle.textContent = isEn ? category.titleEn : category.titleZh;
        catDiv.appendChild(catTitle);

        if (category.noteZh || category.noteEn) {
            var note = document.createElement('p');
            note.className = 'papers-note';
            note.textContent = isEn ? category.noteEn : category.noteZh;
            catDiv.appendChild(note);
        }

        var globalIndex = 1;

        category.periods.forEach(function(period) {
            if (period.titleZh || period.titleEn) {
                var periodTitle = document.createElement('h4');
                periodTitle.className = 'papers-period-title';
                periodTitle.textContent = isEn ? period.titleEn : period.titleZh;
                catDiv.appendChild(periodTitle);
            }

            var ol = document.createElement('ol');
            ol.className = 'papers-list';
            ol.setAttribute('start', globalIndex);

            period.items.forEach(function(text) {
                var li = document.createElement('li');
                li.className = 'paper-item';
                li.innerHTML = formatPaperText(text);
                ol.appendChild(li);
                globalIndex++;
            });

            catDiv.appendChild(ol);
        });

        container.appendChild(catDiv);
    });
}

/* ===== Projects Rendering ===== */

function renderProjects() {
    var container = document.getElementById('projects-container');
    if (!container || typeof projectsData === 'undefined') return;

    var isEn = document.documentElement.lang === 'en';
    container.innerHTML = '';

    projectsData.forEach(function(role) {
        var roleDiv = document.createElement('div');
        roleDiv.className = 'projects-role';

        var roleTitle = document.createElement('h3');
        roleTitle.className = 'projects-role-title';
        roleTitle.textContent = isEn ? role.titleEn : role.titleZh;
        roleDiv.appendChild(roleTitle);

        var ol = document.createElement('ol');
        ol.className = 'projects-list';

        role.items.forEach(function(proj) {
            var li = document.createElement('li');
            li.className = 'project-item';

            var nameSpan = document.createElement('span');
            nameSpan.className = 'project-name';
            nameSpan.textContent = proj.name;
            li.appendChild(nameSpan);

            var meta = document.createElement('span');
            meta.className = 'project-meta';

            var parts = [];
            parts.push(proj.source);
            parts.push(isEn ? 'Period: ' + proj.period : '期間：' + proj.period);
            if (proj.amount) parts.push(proj.amount);

            var grantStr = proj.grant;
            if (proj.internal) grantStr += ' / ' + proj.internal;

            meta.innerHTML = parts.join(' · ') +
                '<br><span class="project-grant">' + grantStr + '</span>';

            li.appendChild(meta);
            ol.appendChild(li);
        });

        roleDiv.appendChild(ol);
        container.appendChild(roleDiv);
    });
}

/* ===== Article Text Expand/Collapse ===== */

function toggleArticleText(el) {
    if (el.classList.contains('collapsed')) {
        el.classList.remove('collapsed');
        el.classList.add('expanded');
    } else {
        el.classList.remove('expanded');
        el.classList.add('collapsed');
    }
}

/* ===== Article PDF Viewer ===== */

function showArticlePdf(filename) {
    var list = document.getElementById('articles-list');
    var viewer = document.getElementById('articles-pdf-viewer');
    if (!list || !viewer) return;
    list.style.display = 'none';
    viewer.style.display = 'block';
    document.getElementById('articles-pdf-frame').src = filename;
    var section = document.getElementById('articles');
    if (section) window.scrollTo({ top: section.offsetTop - 60, behavior: 'smooth' });
}

function hideArticlePdf() {
    var list = document.getElementById('articles-list');
    var viewer = document.getElementById('articles-pdf-viewer');
    if (!list || !viewer) return;
    viewer.style.display = 'none';
    list.style.display = 'block';
    var frame = document.getElementById('articles-pdf-frame');
    if (frame) frame.src = 'about:blank';
}

/* ===== Image Modal ===== */

function showImageModal(src) {
    var overlay = document.getElementById('image-modal-overlay');
    var img = document.getElementById('image-modal-img');
    if (!overlay || !img) return;
    img.src = src;
    overlay.style.display = 'flex';
}

function hideImageModal() {
    var overlay = document.getElementById('image-modal-overlay');
    var img = document.getElementById('image-modal-img');
    if (!overlay || !img) return;
    img.src = '';
    overlay.style.display = 'none';
}

/* ===== Video Modal ===== */

function showVideoModal(youtubeId) {
    if (!/^[a-zA-Z0-9_-]{11}$/.test(youtubeId)) return;
    var overlay = document.getElementById('video-modal-overlay');
    var frame = document.getElementById('video-modal-frame');
    if (!overlay || !frame) return;
    frame.src = 'https://www.youtube.com/embed/' + youtubeId + '?autoplay=1';
    overlay.style.display = 'flex';
}

function hideVideoModal() {
    var overlay = document.getElementById('video-modal-overlay');
    var frame = document.getElementById('video-modal-frame');
    if (!overlay || !frame) return;
    frame.src = 'about:blank';
    overlay.style.display = 'none';
}

/* ===== Initialization ===== */

function initShared() {
    const langBtn = document.querySelector('.language-toggle');
    if (langBtn) langBtn.addEventListener('click', toggleLanguage);

    renderStats();
    renderHonors();
    renderAlumniArticle();
    renderHonorsArticle();
    renderBio();
    renderPapers();
    renderProjects();
    renderNav();
    renderFooter();
    renderBackToTop();
    updateAriaLabels();
}

function renderBackToTop() {
    if (document.querySelector('.back-to-top-shared')) return;

    const isEn = document.documentElement.lang === 'en';
    const btn = document.createElement('button');
    btn.className = 'back-to-top-shared';
    btn.setAttribute('aria-label', isEn ? 'Back to top' : '回到頂部');
    btn.textContent = '↑';

    const main = document.querySelector('main') || document.body;
    main.appendChild(btn);

    window.addEventListener('scroll', function () {
        if (window.pageYOffset > 300) {
            btn.classList.add('visible');
        } else {
            btn.classList.remove('visible');
        }
    }, { passive: true });

    btn.addEventListener('click', function () {
        window.scrollTo({ top: 0, behavior: 'smooth' });
    });
}

/* ===== Microsoft Clarity ===== */

(function(c,l,a,r,i,t,y){
    c[a]=c[a]||function(){(c[a].q=c[a].q||[]).push(arguments)};
    t=l.createElement(r);t.async=1;t.src="https://www.clarity.ms/tag/"+i+"?ref=bwt";
    y=l.getElementsByTagName(r)[0];y.parentNode.insertBefore(t,y);
})(window, document, "clarity", "script", "rzlnthqbys");
