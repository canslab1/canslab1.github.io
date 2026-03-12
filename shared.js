/* ===== Honors Data ===== */

const honorsZh = [
    '<span class="honors-highlight">中國科技大學</span><br>2024.11 榮獲傑出校友<br>2000.11 榮獲優秀校友<br>校史唯一先後獲頒傑出與優秀校友',
    '<span class="honors-highlight">臺北市杏壇芬芳獎</span><br>2024.11 榮獲臺北市教育局公開表揚',
    '<span class="honors-highlight">臺北市老松國民小學學生家長會</span><br>2024.10 榮膺第 34 屆榮譽會長<br>2023.10 擔任臺北市小聯會常務理事<br>2022.10 擔任第 32 和 33 屆會長<br>113 學年起提供畢業生榮譽會長獎，鼓勵多元學習表現優異學童',
    '<span class="honors-highlight">所羅門股份有限公司</span><br>2022.06 榮膺第十二屆董事會獨立董事<br>2022.06 擔任薪酬、審計、永續發展委員',
    '<span class="honors-highlight">中華民國斐陶斐榮譽學會獎</span><br>2006.06 榮獲國立交通大學公開頒發<br>表彰卓越學術成就與優異品德操守',
    '<span class="honors-highlight">臺北市杏壇芬芳獎推薦文</span><br>源源不絕無私大愛，守護美好優質老松'
];

const honorsEn = [
    '<span class="honors-highlight">China University of Technology</span><br>2024.11 Distinguished Alumni Award<br>2000.11 Outstanding Alumni Award<br>the only recipient in university history of both Distinguished and Outstanding Alumni honors.',
    '<span class="honors-highlight">2024.10 Taipei Honored Contributors to Education</span><br>Bottomless Devotion to Taipei Lao-Song Elementary School',
    '<span class="honors-highlight">Laosong Elementary School, Taipei City</span><br>2024.10 Honorary President, 34th Parents&#39; Association<br>2022.10-2024.10 President of the 32nd and 33rd Parents&#39; Association<br>Starting from the 113th Academic Year, the Honorary PTA President\'s Award will be presented to graduating students',
    '<span class="honors-highlight">Solomon Technology Corporation</span><br>2022.06 Served as Independent Director on the 12th Board of Directors of Solomon Co., Ltd.<br>2022.06 Served on the Compensation, Audit, and Sustainability Committees of the 12th Board of Directors of Solomon Co., Ltd.',
    '<span class="honors-highlight">Phi Tau Phi Honor Award</span><br>2006.06 Awarded by National Chiao Tung University,<br>in recognition of academic excellence and exemplary character',
    '<span class="honors-highlight">Recommendation Article for the Taipei Honored contributors to education</span>'
];

/* ===== Navigation Data ===== */

const navItemsZh = [
    { label: '介紹', href: null, section: 'overview' },
    { label: '實驗室', href: null, section: 'lab' },
    { label: '程式', href: null, section: 'software' },
    { label: '著作', href: 'https://scholar.google.com/citations?user=0klfzfAAAAAJ&hl=en', embed: true },
    { label: '計畫', href: 'https://arspb.nstc.gov.tw/NSCWebFront/modules/talentSearch/talentSearch.do?action=initRsm17new&rsNo=2a2ef24adef742f58349c3c533bb7402&LANG=chi', embed: true },
    { label: '履歷', href: 'https://canslab1.github.io/CV.pdf', embed: true },
    { label: '維基', href: 'https://sites.google.com/view/gscott-huang' },
    { label: '臉書', href: 'https://www.facebook.com/gscott.huang/' },
    { label: '故事', href: 'https://canslab1.github.io/stories.html' },
    { label: '相簿', href: 'https://www.dropbox.com/scl/fo/ejwke56gkn1erv7meid4k/AHCPQXAw_RpiEn_PeiBfCuI?rlkey=b8k879fufdh98x2s3fzkatk5y&e=1&dl=0' },
    { label: '評論', href: 'https://www.google.com/search?q=%22Huang+Chung-yuan%22+%22Taipei+Times%22+%22Chang+Gung%22+site%3Awww.taipeitimes.com' },
    { label: '投書', href: 'https://www.google.com/search?q=%22%E9%BB%83%E5%B4%87%E6%BA%90%22+%22%E8%87%AA%E7%94%B1%E6%99%82%E5%A0%B1%22+site%3Atalk.ltn.com.tw' }
];

const navItemsEn = [
    { label: 'About', href: null, section: 'overview' },
    { label: 'Lab', href: null, section: 'lab' },
    { label: 'Software', href: null, section: 'software' },
    { label: 'Papers', href: 'https://scholar.google.com/citations?user=0klfzfAAAAAJ&hl=en', embed: true },
    { label: 'Projects', href: 'https://arspb.nstc.gov.tw/NSCWebFront/modules/talentSearch/talentSearch.do?action=initRsm17new&rsNo=2a2ef24adef742f58349c3c533bb7402&LANG=chi', embed: true },
    { label: 'CV', href: 'https://canslab1.github.io/CV.pdf', embed: true },
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

    renderStats();
    renderHonors();
    renderHonorsArticle();
    renderNav();
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

/* ===== Render Functions ===== */

function renderStats() {
    const container = document.getElementById('stats-container');
    if (!container) return;
    const isEn = document.documentElement.lang === 'en';
    const stats = isEn ? statsEn : statsZh;
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

var _activeSection = 'overview';

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

    items.forEach(item => {
        if (item.section) {
            const btn = document.createElement('button');
            btn.className = 'nav-item' + (item.section === _activeSection ? ' active' : '');
            btn.textContent = item.label;
            btn.addEventListener('click', (e) => showSection(item.section, e));
            container.appendChild(btn);
        } else if (item.embed) {
            const btn = document.createElement('button');
            btn.className = 'nav-item' + (_activeSection === 'embed' ? '' : '');
            btn.textContent = item.label;
            btn.addEventListener('click', (e) => showEmbed(item.href, e));
            container.appendChild(btn);
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

/* ===== Honors Article Data ===== */

const honorsArticleZh = [
    '「蘇○○，早安！何○○，早安！」',
    '如冬日暖陽般親切地喊出學童姓名的崇源會長，每天早晨總是站在穿堂旁學生家長會辦公室門口，以充滿朝氣的笑容，和藹可親的問候，迎接每一位踏入老松校園的學生。在學童的眼中，他是永遠面帶微笑，樂意彎下腰、蹲下身，傾聽他們說話，關懷他們心情的學生家長會長。',
    '崇源會長是老松校友，畢業於民國七十一年，家族三代都就讀老松，舅舅更長期擔任老松國小訓導主任。崇源會長常掛在嘴邊的一句話：「老松國小的大小事，都是我的大事，能夠協助母校親師生，是畢生的榮幸！」',
    '崇源會長帶著一顆感恩、奉獻、回饋的心，投入老松學生家長會志工服務的行列，期間經歷了家長會委員、志工總幹事、活動組長、常務委員、副會長等職務，非常熟悉會務。在疫情接近尾聲之際，崇源會長結合了一群來自低、中、高年級，橫跨園藝、保健、導護、愛閱、圖書、ＥＱ、學輔所有組別的熱忱家長志工，秉持著自發、互動、共好的教育理念，入校服務，全力支援老松國小於疫後陸續恢復的各項活動，積極協助學童在花木扶疏、綠意盎然的校園環境中探索與學習。',
    '時間是最寶貴的資源。為了協助疫後重返校園的老松師生與愛心志工，崇源會長毫無私心，奉獻出寶貴的時間與無比的愛心，每日坐鎮學生家長會辦公室，用充滿活力與正向積極的任事態度，激勵團隊夥伴，透過集思廣益，淬鍊出多元的趣味創意與靈活想法，盡心盡力辦好學生家長會務，協助校務與教育工作的發展，關心老松孩童的校園生活。',
    '曾是國小校警，如今擔任保全的蕭○○先生，常對崇源會長豎起大姆指說：「會長，您真的很厲害！整日守護在學校，處理會務，關懷親師生，其他學校的家長會辦公室通常大門深鎖很少使用，只有開會和活動才有人氣，兩者差異很大。」',
    '除了獻出寶貴的時間，崇源會長也以大學教授平日從事研究的嚴謹態度，檢視學生家長會辦理的親子成長、親職講座、志工成長、感恩餐會、一日遊等活動。從前置準備，到執行細節，再到預期結果，他總是不厭其煩和團隊夥伴反覆推敲活動的每項細節，可能的意外，以及臨場可行的補救之道。崇源會長常對團隊夥伴說：「舉辦活動，就要辦得熱鬧有意義，讓參與活動的親師生與愛心志工收穫滿滿，身心靈滿載而歸。」',
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

const honorsFigures = [
    {
        src: 'https://canslab1.github.io/images/honors1.png',
        alt: '老松國小．學生家長會．黃崇源會長過往二年的部分成果展現',
        captionZh: '▲ 在老松學童的眼中，崇源會長就是聖誕老公公，帶來歡笑，帶來令人難忘的回憶。',
        captionEn: '▲ To the children of Lao-Song Elementary, President Chung-Yuan is their Santa Claus\u2014bringing joy, spreading laughter, and leaving behind memories that linger long after the moment has passed.'
    },
    {
        src: 'https://canslab1.github.io/images/honors2.png',
        alt: '老松國小．學生家長會．黃崇源會長過往二年的部分成果展現',
        captionZh: '▲ 崇源會長出錢出力，盡心盡力辦好學生家長會主辦的親師生及愛心志工活動。',
        captionEn: '▲ President Chung-Yuan contributed both time and resources with wholehearted dedication, ensuring that every event for teachers, students, and volunteers was carried out with care, creativity, and purpose.'
    },
    {
        src: 'https://canslab1.github.io/images/honors3.png',
        alt: '老松國小．學生家長會．黃崇源會長過往二年的部分成果展現',
        captionZh: '▲ 在崇源會長付出寶貴的時間與無比的愛心後，老松國小親師生同心協力共好共榮。',
        captionEn: '▲ Through the generous time and boundless care of President Chung-Yuan, the Lao-Song community\u2014teachers, parents, and students alike\u2014came together in harmony, working hand in hand toward shared growth and collective flourishing.'
    },
    {
        src: 'https://canslab1.github.io/images/honors4.png',
        alt: '老松國小．學生家長會．黃崇源會長過往二年的部分成果展現',
        captionZh: '▲ 崇源會長永遠面帶微笑，和老松所有學童（包含幼兒園幼童）相處融洽。',
        captionEn: '▲ With a warm and ever-present smile, President Chung-Yuan connects effortlessly with every child at Lao-Song\u2014from kindergarteners to upper-grade students\u2014creating bonds built on kindness, trust, and joy.'
    }
];

function renderHonorsArticle() {
    const container = document.getElementById('honors-article');
    if (!container) return;

    const isEn = document.documentElement.lang === 'en';
    const paragraphs = isEn ? honorsArticleEn : honorsArticleZh;

    container.innerHTML = '';

    const p = document.createElement('p');
    p.className = 'honors-paragraph';
    paragraphs.forEach((text, i) => {
        p.appendChild(document.createElement('br'));
        p.appendChild(document.createTextNode('\u2003\u2003' + text));
        if (i < paragraphs.length - 1) p.appendChild(document.createElement('br'));
    });
    container.appendChild(p);

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

/* ===== Initialization ===== */

function initShared() {
    const langBtn = document.querySelector('.language-toggle');
    if (langBtn) langBtn.addEventListener('click', toggleLanguage);

    renderStats();
    renderHonors();
    renderHonorsArticle();
    renderNav();
    renderFooter();
}

/* ===== Microsoft Clarity ===== */

(function(c,l,a,r,i,t,y){
    c[a]=c[a]||function(){(c[a].q=c[a].q||[]).push(arguments)};
    t=l.createElement(r);t.async=1;t.src="https://www.clarity.ms/tag/"+i+"?ref=bwt";
    y=l.getElementsByTagName(r)[0];y.parentNode.insertBefore(t,y);
})(window, document, "clarity", "script", "rzlnthqbys");
