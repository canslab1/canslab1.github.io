# canslab1.github.io

黃崇源教授（Prof. Chung-Yuan Huang, Ph.D.）個人學術網站，託管於 GitHub Pages。

> **最後更新**：2026-04-02 ｜ **更新歷史**：[Commit Log](https://github.com/canslab1/canslab1.github.io/commits/master)

## 目錄

- [網站概覽](#網站概覽)
- [快速開始（本機開發）](#快速開始本機開發)
- [教授簡歷](#教授簡歷)
- [頁面結構](#頁面結構)
- [首頁（index.html）內容區段](#首頁indexhtml內容區段)
- [著作列表（papers-data.js）](#著作列表papers-datajs)
- [計畫列表（projects-data.js）](#計畫列表projects-datajs)
- [研究計畫摘要（CV.pdf 內容）](#研究計畫摘要cvpdf-內容)
- [獎項與榮譽完整列表](#獎項與榮譽完整列表)
- [技術架構](#技術架構)
- [部署與發佈流程](#部署與發佈流程)
- [瀏覽器相容性](#瀏覽器相容性)
- [家族故事頁面（stories.html）](#家族故事頁面storieshtml)
- [404 錯誤頁面](#404-錯誤頁面)
- [SEO 設定](#seo-設定)
- [LLM 搜尋設定](#llm-搜尋設定)
- [IndexNow 即時索引](#indexnow-即時索引)
- [子 Repo 架構與連結關係](#子-repo-架構與連結關係)
- [分析工具與搜尋引擎管理](#分析工具與搜尋引擎管理)
- [第三方外部資源](#第三方外部資源)
- [無障礙設計（Accessibility）](#無障礙設計accessibility)
- [網站地圖（sitemap.xml）](#網站地圖sitemapxml)
- [外部連結](#外部連結)
- [內容管理指南](#內容管理指南)
- [常見維護操作](#常見維護操作)
- [授權說明](#授權說明)
- [貢獻指南](#貢獻指南)
- [備份與災難復原](#備份與災難復原)
- [疑難排解（Troubleshooting）](#疑難排解troubleshooting)

## 網站概覽

| 項目 | 說明 |
|------|------|
| 網址 | https://canslab1.github.io/ |
| 擁有者 | 黃崇源教授（長庚大學資訊工程學系） |
| 語言 | 中文 / English（即時切換） |
| 託管 | GitHub Pages |
| 電郵 | gscott@mail.cgu.edu.tw |

## 快速開始（本機開發）

本網站為純靜態網站，無需安裝任何套件。

```bash
# 1. 取得原始碼
git clone https://github.com/canslab1/canslab1.github.io.git
cd canslab1.github.io

# 2. 啟動本機伺服器（任選一種）
python3 -m http.server 8000           # Python
npx serve .                           # Node.js（需先安裝 Node.js）
open index.html                       # 或直接用瀏覽器開啟（部分功能可能受限於 CORS）

# 3. 瀏覽
# 前往 http://localhost:8000
```

> **注意**：直接用瀏覽器開啟 `index.html`（`file://` 協議）時，某些瀏覽器可能因 CORS 限制導致 JS 載入失敗，建議使用本機伺服器。

## 教授簡歷

### 現任職務
- **長庚大學** 工學院 資訊工程學系 教授（教字第 141150 號，104.08 起）
- **長庚大學** 智慧運算學院 / 人工智慧學系 合聘教授（109.10 起）
- **老松國小** 115 年校務顧問（115.03 起）
- **老松國小** 學生家長會 榮譽會長（113.10 起）

### 研究領域
- 複雜網絡電腦建模與模擬
- 社會現象電腦建模與模擬
- 傳染病傳播動態電腦模擬

### 主持實驗室
複雜適應性網絡與系統實驗室（Complex Adaptive Networks and Systems Lab, CANS Lab）

### 教育背景

| 學位 | 學校 | 系所 | 期間 |
|------|------|------|------|
| 博士 | 交通大學 | 資訊科學研究所 | 89.08–94.11 |
| 碩士 | 交通大學 | 資訊科學研究所 | 87.08–89.07 |
| 學士 | 中國科技大學 | 電子資料處理科 | 74.08–81.07 |

### 教職經歷

| 職稱 | 單位 | 期間 |
|------|------|------|
| 教授 | 長庚大學 資訊工程學系 | 104.08 至今 |
| 副教授 | 長庚大學 資訊工程學系 | 99.08–104.07 |
| 助理教授 | 長庚大學 資訊工程學系 | 95.08–99.07 |
| 助理教授兼系主任 | 元培科技大學 資訊工程學系 | 94.08–95.07 |
| 專任講師 | 中國科技大學 資訊工程學系 | 91.08–94.07 |

### 業界經歷（攻讀博士前，共 11 年）
- 宇博電腦（新加坡）— 資深系統設計師
- 環世股份有限公司 — 系統工程師
- 百點連鎖資訊教室 — 教學部主任 / 總經理特別助理
- 鼎眾電腦股份有限公司 — 軟體開發部主任

### 授課科目
- 大學部：網路應用軟體設計（IT3019）、物件導向軟體設計（IT3004）、系統程式（IT2017）
- 研究所：機器學習（ITM114）、資料探勘（ITM063）、柔性計算（ITM141）

## 頁面結構

本網站採用**單頁式架構**（Single Page Application），`index.html` 包含八個內容區段（介紹 / 榮譽 / 綜覽 / 著作 / 計畫 / 程式 / 評論）以及內嵌 PDF 檢視（履歷），透過導覽列按鈕切換，不重新載入頁面。

| 檔案 | 頁面 | 說明 |
|------|------|------|
| `index.html` | 首頁（單頁式） | 包含八個區段：介紹（簡短自傳）、榮譽（責任與肯定 + 傑出校友推薦文 + 杏壇芬芳錄推薦文）、綜覽（28 張學術數據卡片）、著作（著作目錄）、計畫（計畫列表）、程式（開源程式）、評論（媒體投書），以及內嵌 PDF 檢視（履歷） |
| `stories.html` | 家族故事 | 個人家族故事集 |
| `404.html` | 錯誤頁 | 自訂 404 頁面 |
| `CYHuang.html` | 重導向頁 | 已移至 `index.html#papers`（`noindex`、`meta refresh`、`canonical`） |
| `software.html` | 重導向頁 | 已移至 `index.html#software`（`noindex`、`meta refresh`、`canonical`） |
| `lab.html` | 重導向頁 | 已移至 `index.html#lab`（`noindex`、`meta refresh`、`canonical`） |
| `CV.pdf` | 履歷 | 完整學術履歷（PDF） |
| `shared.css` | 共用樣式 | 各頁面共用的 CSS 樣式表 |
| `shared.js` | 共用腳本 | 渲染函式、語言切換、導覽列/Footer、區段切換 |
| `stats-data.js` | 卡片資料 | 28 張學術數據卡片的中英文資料 |
| `papers-data.js` | 論文資料 | 171 篇著作列表（5 大類，中英文） |
| `projects-data.js` | 計畫資料 | 37 件研究計畫列表（3 種角色，中英文） |
| `build-prerender.py` | 建構工具 | 解析 JS 資料檔，產生預渲染 HTML 及 JSON-LD 結構化資料 |
| `css/project-page.css` | 子 repo 共用樣式 | 8 個子 repo 專案頁面共用 CSS（header、markdown-body、footer、響應式） |
| `js/readme-loader.js` | 子 repo 共用腳本 | 自動偵測 repo 名稱、fetch README.md、渲染 Markdown 與語法高亮 |
| `stories.css` | 故事樣式 | 家族故事頁面專屬 CSS |
| `stories.js` | 故事腳本 | 家族故事頁面專屬 JS |

## 首頁（index.html）內容區段

### Header
- 大頭照（`images/IMG-2.jpg`）
- 姓名：長庚大學 黃崇源 教授
- 副標：複雜適應性網絡與系統實驗室 主持人

### 導覽列

中文：介紹 | 榮譽 | 綜覽 | 著作 | 計畫 | 程式 | 評論 | 履歷 | 故事 | 相簿

英文：About | Honors | Overview | Papers | Projects | Softwares | Op-Eds | CV | Stories | Photos

前七項（介紹 / 榮譽 / 綜覽 / 著作 / 計畫 / 程式 / 評論）為頁面內區段切換，履歷為嵌入式 PDF 檢視，其餘兩項（故事 / 相簿）為外部連結。外部連結使用具名 `target`（如 `canslab-故事`），重複點擊同一連結時瀏覽器會跳到已開啟的視窗而非再開新分頁。

### 區段一：簡短自傳（介紹，預設顯示）

顯示黃崇源教授的中英文簡短自傳，內容涵蓋教育背景、教職經歷、研究興趣及學術成就。

> **維護提示**：自傳內容更新請修改 `shared.js` 中的 `bioZh[]` 和 `bioEn[]` 陣列。

### 區段二：責任與肯定（榮譽）

包含三個部分：

**一、責任與肯定（榮譽列表，共 6 大項）**
- **中國科技大學** — 2024 傑出校友（校史唯一先後獲頒傑出與優秀校友）— 附「專訪」影片按鈕（YouTube 彈窗）、「照片一」「照片二」按鈕
- **臺北市杏壇芬芳獎** — 2024 榮獲臺北市教育局公開表揚 — 附「照片一」「照片二」按鈕 — 源源不絕無私大愛，守護美好優質老松
- **臺北市老松國民小學** — 2026 榮膺 115 年校務顧問（附照片）、2025 教師節捐款伍拾萬元整（附照片一至五）、111 學年起連續四學年捐款總額超過壹佰萬元整
- **臺北市老松國民小學學生家長會** — 2024 第 34 屆榮譽會長（附照片）、2023 臺北市小聯會常務理事（附照片一至四）、2022 第 32 和 33 屆會長、113 學年起設立畢業生榮譽會長獎（附照片）
- **所羅門股份有限公司** — 2022 獨立董事（薪酬、審計、永續發展委員）— 附「照片一」「照片二」按鈕
- **中華民國斐陶斐榮譽學會獎** — 2006（國立交通大學）— 附「照片一」「照片二」按鈕

**二、老松國小傑出校友推薦文**
收錄完整的臺北市老松國小 130 週年校慶傑出校友貢獻獎推薦文（中英文）。

**三、杏壇芬芳錄推薦文**
收錄完整的臺北市杏壇芬芳獎推薦文（中英文），記述崇源會長於老松國小學生家長會的奉獻，以及四張活動照片（`honors1–4.png`）。

> **維護提示**：榮譽列表更新請修改 `shared.js` 中的 `honorsZh[]` 和 `honorsEn[]` 陣列；傑出校友推薦文更新請修改 `alumniArticleZh[]` 和 `alumniArticleEn[]` 陣列；杏壇芬芳錄推薦文更新請修改 `honorsArticleZh[]` 和 `honorsArticleEn[]` 陣列。

### 區段三：學術數據卡片（綜覽）

直接顯示 28 張數據卡片（無標題、無綠線分隔），涵蓋：
- 期刊論文：44 篇
- 國際研討會：50 篇
- 專書專章：10 篇
- 國內研討會：19 篇
- 文章及採訪：49 次
- 計畫主持人：22 次（國科會預算 1,408 萬）
- 共同主持人：14 次
- Google Scholar：Citations 994、h-index 18、i10-index 30
- 指導學生：專題生 51 人、碩士 25 人、博士 2 人
- 學術服務：論文審稿 68 篇、會議審稿 27 次、議程委員 77 次、計劃審查 11 次、學術審查 10 次、專題演講 35 次
- 教職年資：教授 10 年、副教授 5 年、助理教授 5 年、講師 5 年
- 其他經歷：業界 11 年、合聘 4 單位、獨立董事 3 年、學生家長會 8 年

> **維護提示**：數據更新請修改 `stats-data.js` 中的 `statsZh[]` 和 `statsEn[]` 陣列，不需異動 `index.html` 或 `shared.js`。

### 區段四：著作目錄（著作）

顯示黃崇源教授 171 篇著作，分為 5 大類別，詳見下方「著作列表」章節。

> **維護提示**：新增或修改論文只需編輯 `papers-data.js`，然後執行 `python3 build-prerender.py` 重新產生預渲染 HTML。

### 區段五：計畫列表（計畫）

顯示黃崇源教授 37 件研究計畫，分為 3 種角色，詳見下方「計畫列表」章節。

> **維護提示**：新增或修改計畫只需編輯 `projects-data.js`，然後執行 `python3 build-prerender.py` 重新產生預渲染 HTML。

### 區段六：開源程式（程式）

展示 CANS 實驗室開發的 8 套開源研究軟體，每套皆附論文縮圖、應用截圖與 GitHub 連結。所有程式碼以 MIT 授權公開。

| 軟體 | 說明 | 語言 | GitHub | README |
|------|------|------|--------|--------|
| **EpiRank** | 基於非對稱通勤網絡的疫情風險分析 | Python | [canslab1/EpiRank](https://github.com/canslab1/EpiRank) | [README](https://canslab1.github.io/EpiRank/) |
| **MV17** | 以 K-core 為基礎的多屬性節點重要性排序 | Python | [canslab1/MV17](https://github.com/canslab1/MV17) | [README](https://canslab1.github.io/MV17/) |
| **CASMIM** | 結合細胞自動機與社會鏡像身份的 SARS 疫情模擬 | Python | [canslab1/CASMIM](https://github.com/canslab1/CASMIM) | [README](https://canslab1.github.io/CASMIM/) |
| **HETA** | 無參數複雜網絡邊分類 | Python | [canslab1/HETA](https://github.com/canslab1/HETA) | [README](https://canslab1.github.io/HETA/) |
| **HATA** | 有向網絡弧分類工具 | Python | [canslab1/HATA](https://github.com/canslab1/HATA) | [README](https://canslab1.github.io/HATA/) |
| **BCAT** | 有界信心意見動態與採用門檻創新擴散混合模擬模型 | Python / NetLogo | [canslab1/BCAT](https://github.com/canslab1/BCAT) | [README](https://canslab1.github.io/BCAT/) |
| **SRAC-Agent** | 演化空間囚徒困境賽局中自我聲譽感知機制模擬器 | Python | [canslab1/SRAC-Agent](https://github.com/canslab1/SRAC-Agent) | [README](https://canslab1.github.io/SRAC-Agent/) |
| **AED2** | 使用基因演算法最佳化 AED 配置地點 | C++ / Python | [canslab1/AED2](https://github.com/canslab1/AED2) | [README](https://canslab1.github.io/AED2/) |

> **維護提示**：軟體項目更新請直接修改 `index.html` 中「程式」區段（`<section id="software">`）的 `.software-card` 區塊。

### 區段七：媒體投書（評論）

收錄黃崇源教授發表於國內外報紙的投書文章共 13 篇，按報社分組（Taipei Times 5 篇、自由時報 5 篇、風傳媒 3 篇），同一報社內按日期由新到舊排列。每篇文章顯示標題（含外部連結）、發表日期與作者。

> **維護提示**：新增投書請直接修改 `index.html` 中「評論」區段（`<section id="press">`）的 `.software-card` 區塊，按報社分組並維持日期由新到舊排序。

### 回到頂部按鈕

當使用者向下捲動超過 300px 時，頁面右下角會出現一個綠底白色向上箭頭的圓形按鈕（48px），點擊後平滑捲動回頁面頂部。此按鈕定位在距離底部 75px 處，避免被 Footer 遮擋。

> **維護提示**：按鈕樣式在 `shared.css` 的 `.back-to-top-shared`，邏輯在 `shared.js` 的 `renderBackToTop()` 函式中。`stories.html` 另有獨立實作（紫底按鈕）。

### Footer

頁面底部固定列，顯示 9 個 logo 圖示連結。滑鼠移到任一 logo 上時，logo 會放大並顯示綠色提示框，標示即將連結到的網站名稱。

| # | Logo | 連結目標 |
|---|------|---------|
| 1 | 長庚資工系 | cgu.edu.tw/csie |
| 2 | 老松國小 | tlsps.tp.edu.tw |
| 3 | 杏壇芬芳錄 | xingtan.tiec.tp.edu.tw |
| 4 | ORCiD | orcid.org |
| 5 | Facebook | facebook.com/gscott.huang |
| 6 | Google Sites | sites.google.com/view/gscott-huang |
| 7 | Google Scholar | scholar.google.com |
| 8 | GitHub | github.com/canslab1 |
| 9 | 長庚 Pure | pure.lib.cgu.edu.tw |

Footer 使用縮圖（`*-thumb.*`）以減少頁面載入大小。Hover 時顯示自訂綠色漸層提示框（`.footer-tooltip`），帶向下箭頭指向 logo。

> **維護提示**：Footer 連結更新請修改 `shared.js` 中的 `footerLinks[]` 陣列。如新增 logo 請同時製作縮圖版本（建議寬度 128px 以內）。提示框文字自動取自 `label` 欄位並去除「前往」前綴。

## 著作列表（papers-data.js）

首頁「著作」區段顯示黃崇源教授 171 篇著作，資料儲存於 `papers-data.js`，分為 5 大類別：

| # | 類別 | 篇數 | 期別分組 |
|---|------|------|---------|
| 1 | SCI / SSCI / EI 等級期刊論文 | 44 | 依教職時期分 5 組 |
| 2 | 國際會議論文 | 49 | 單一列表 |
| 3 | 專章、學位及其它論文 | 10 | 單一列表 |
| 4 | 國內外媒體投書及採訪 | 49 | 單一列表 |
| 5 | 國內研討會論文 | 19 | 單一列表 |

- 論文內容會由 `build-prerender.py` 預渲染為靜態 HTML，嵌入 `index.html` 的 `<div id="papers-container">`
- 瀏覽器載入後 `shared.js` 的 `renderPapers()` 會重新渲染（支援語言切換）
- 作者名 "Huang C.Y." 自動加粗，帶 DOI 的論文自動產生超連結

## 計畫列表（projects-data.js）

首頁「計畫」區段顯示黃崇源教授 37 件研究計畫，資料儲存於 `projects-data.js`，分為 3 種角色：

| 角色 | 件數 | 包含資訊 |
|------|------|---------|
| 主持人 | 22 | 期間、來源、計畫編號、校內編號、經費、名稱 |
| 共同主持人 | 14 | 期間、來源、計畫編號、名稱 |
| 研究人員 | 1 | 期間、來源、計畫編號、名稱 |

- 計畫內容由 `build-prerender.py` 預渲染為靜態 HTML，嵌入 `index.html` 的 `<div id="projects-container">`
- 瀏覽器載入後 `shared.js` 的 `renderProjects()` 會重新渲染（支援語言切換）

## 研究計畫摘要（CV.pdf 內容）

### 主持人計畫（22 件，國科會/科技部）

| 期間 | 計畫名稱 | 預算 |
|------|--------|------|
| 109–110 | 以地理資訊系統為基礎之人工智慧基因演算法暨時空權重雙模型 | 113 萬 |
| 107–108 | 整合資訊熵值及局域拓樸特性來識別複雜網絡超級傳播者 | 78.4 萬 |
| 106–107 | 使用規則式拓樸策略來快速識別複雜網絡社群及其階層結構 | 53.9 萬 |
| 105–106 | 線上熟人社群網絡演化之交友資源與記憶的影響 | 61.4 萬 |
| 104–105 | 反覆式囚犯困局中可自我聲譽調整能力之智慧型代理人的影響 | 82.1 萬 |
| 103–104 | 最適新型流感交通阻絕策略之基因演算法優選與防疫成本效益分析 | 50 萬 |
| 100–101 | 台灣新型流感 H1N1 疫情之計量地理分析與時空傳播網絡監測平台 | 50.7 萬 |
| 99–100 | 建構「台灣通勤人口暨交通運輸網絡模型」 | 78 萬 |
| 98–99 | 居家隔離檢疫政策成本效益及最佳實施策略評估分析 | 90 萬 |
| 95–96 | 小世界流行病學建模與公共衛生政策之評估 | 44.9 萬 |

### 共同主持人計畫（14 件）
涵蓋傳染病監測、藥物濫用預警、GPU 演算法、P2P 網路等領域。

### 合作機構
國科會（NSTC/MOST）、新興病毒感染研究中心（CMRPD）、疾病管制局（DOH）、食品藥物管理局（TFDA）

## 獎項與榮譽完整列表

- **2026** 臺北市老松國民小學 115 年校務顧問
- **2025** 111 學年起，連續四學年捐款總額超過新臺幣壹佰萬元整（老松國小）
- **2025** 教師節捐款新台幣伍拾萬元整予老松國小
- **2024** 臺北市杏壇芬芳獎（老松國小推薦，含頒獎照片）— 源源不絕無私大愛，守護美好優質老松
- **2024** 中國科技大學 59 週年校慶傑出校友（含專訪影片、活動照片）
- **2019** 第 17 屆台塑企業應用技術研討會研發論文暨海報競賽學校組優勝第一名
- **2016** ICIAE 2016 Best Student Paper Award（指導學生傅昱翔）
- **2015** 長庚大學工學院校外實習課程成果競賽佳作
- **2014** Simulation 期刊 Editor's Choice 排行榜第二名
- **2013** e-CASE & e-Tech 2013 最佳論文首獎 Shieh Tung-Min Memorial Award
- **2012** 中華電信創新應用大賽季軍
- **2006** 中華民國斐陶斐榮譽學會 2006 年度榮譽會員
- **2001** 第二屆台灣工業銀行創業大賽佳作
- **2000** 中國科技大學第 12 屆優秀校友
- 研發全國首套社會網絡式傳染病傳播動態電腦模擬系統 CASMIM

## 技術架構

### 檔案結構
```
canslab1.github.io/
├── index.html           # 首頁（單頁式，含介紹/榮譽/綜覽/著作/計畫/程式/評論七個區段 + 內嵌 PDF）
├── stories.html         # 家族故事
├── stories.css          # 故事頁面專屬 CSS
├── stories.js           # 故事頁面專屬 JS
├── shared.css           # 共用 CSS 樣式（含回到頂部按鈕、Footer 提示框）
├── shared.js            # 共用 JS 渲染函式與區段切換（含回到頂部按鈕邏輯）
├── stats-data.js        # 28 張學術數據卡片資料
├── papers-data.js       # 171 篇著作列表資料
├── projects-data.js     # 37 件研究計畫列表資料
├── build-prerender.py   # 預渲染建構工具（產生靜態 HTML 及 JSON-LD）
├── 404.html             # 自訂錯誤頁
├── CYHuang.html         # 重導向頁（→ index.html#papers）
├── software.html        # 重導向頁（→ index.html#software）
├── lab.html             # 重導向頁（→ index.html#lab）
├── CV.pdf               # 學術履歷
├── d22a81b3...52.txt    # IndexNow API key 驗證檔
├── security.txt         # RFC 9116 安全聯絡資訊
├── manifest.json        # PWA Web App Manifest（雙語描述）
├── feed.xml             # Atom feed（媒體投書 / 文章動態）
├── humans.txt           # 網站製作者資訊
├── README.md            # 本文件
├── LICENSE              # MIT 授權條款
├── css/
│   └── project-page.css # 8 個子 repo 專案頁面共用 CSS
├── js/
│   └── readme-loader.js # 8 個子 repo README 渲染共用 JS
├── .github/workflows/
│   └── indexnow.yml    # IndexNow 自動提交 workflow
├── .gitignore          # Git 忽略規則（.DS_Store、.claude/）
├── robots.txt          # 爬蟲規則
├── sitemap.xml         # 網站地圖
├── llms.txt            # LLM 搜尋摘要
├── BingSiteAuth.xml    # Bing 驗證
├── google9ba7e...html  # Google 驗證
├── ahrefs_b1f...       # Ahrefs 驗證
├── openai-domain-...   # OpenAI 驗證
└── images/
    ├── IMG-2.jpg       # 教授大頭照
    ├── CANS.png        # CANS Lab logo
    ├── Photos1-1.webp  # 故事頁面用照片（已轉換為 WebP）
    ├── Photos1-2.webp  # 故事頁面用照片（已轉換為 WebP）
    ├── csie.png        # 長庚資工系 logo（原圖，用作 favicon）
    ├── csie-thumb.png  # 長庚資工系 logo（Footer 縮圖）
    ├── cgu-thumb.png   # 長庚大學 logo（Footer 縮圖，用於 Pure 連結）
    ├── laosong-thumb.png  # 老松國小 logo（Footer 縮圖）
    ├── favicon-thumb.ico  # 杏壇芬芳錄 logo（Footer 縮圖）
    ├── facebook-thumb.png # Facebook logo（Footer 縮圖）
    ├── google-thumb.png   # Google logo（Footer 縮圖）
    ├── google-scholar-thumb.png # Google Scholar logo（Footer 縮圖）
    ├── github-thumb.png   # GitHub logo（Footer 縮圖）
    ├── honors1–4.webp  # 杏壇芬芳獎照片（已轉換為 WebP）
    ├── honors/         # 榮譽活動照片（彈窗顯示，全部 WebP 格式）
    │   ├── laosong-advisor.webp         # 老松國小 115 年校務顧問聘書
    │   ├── cute-alumni-group.webp       # 中國科大傑出校友大會團體照
    │   ├── cute-distinguished-alumni.webp # 中國科大傑出校友合影
    │   ├── xtff-award.webp              # 杏壇芬芳獎頒獎照
    │   ├── xtff-group.webp              # 杏壇芬芳獎團體合影
    │   ├── laosong-honorary-president.webp # 老松國小榮譽會長運動會照
    │   ├── laosong-president-award.webp # 老松國小榮譽會長獎頒獎照
    │   ├── laosong-donation-1.webp      # 教師節捐款照片一
    │   ├── laosong-donation-2.webp      # 教師節捐款照片二
    │   ├── laosong-donation-3.webp      # 教師節慶祝活動海報
    │   ├── laosong-donation-4.webp      # 捐款收支明細表
    │   ├── laosong-donation-5.webp      # 捐款支票照片
    │   ├── laosong-xlh-delegates.webp   # 小聯會代表大會合照
    │   ├── laosong-xlh-board.webp       # 小聯會理監事會議
    │   ├── laosong-xlh-election.webp    # 小聯會常務理事選舉結果
    │   ├── laosong-xlh-west-district.webp # 小聯會西區理事選舉結果
    │   ├── solomon-newyear.webp         # 所羅門集團新年晚會
    │   ├── solomon-shareholders.webp    # 所羅門股東常會
    │   ├── ptphs-family.webp            # 斐陶斐頒獎與家人合影
    │   └── ptphs-ceremony.webp          # 斐陶斐授證典禮
    └── software/       # 研究軟體截圖與論文縮圖
        ├── epirank-*.png   # EpiRank 相關圖片
        ├── mv17-*.png      # MV17 相關圖片
        ├── casmim-*.png    # CASMIM 相關圖片
        ├── heta-*.png      # HETA 相關圖片
        ├── hata-*.png      # HATA 相關圖片
        ├── bcat-*.png      # BCAT 相關圖片
        ├── srac-*.png      # SRAC-Agent 相關圖片
        └── aed2-*.png      # AED2 相關圖片
```

### 共用檔案說明

#### `css/project-page.css`
8 個子 repo 專案頁面共用 CSS，詳見「[子 Repo 架構與連結關係](#子-repo-架構與連結關係)」章節。

#### `js/readme-loader.js`
8 個子 repo README 渲染共用 JS，詳見「[子 Repo 架構與連結關係](#子-repo-架構與連結關係)」章節。

#### `shared.css`
全站共用 CSS，包含：
- Reset 樣式、Layout 容器
- Header、導覽列、Footer 樣式
- Footer 提示框（`.footer-tooltip`）— 綠色漸層背景，hover 時淡入顯示
- 學術數據卡片（`.stat-card`）、榮譽列表（`.honors-list`）
- 軟體卡片（`.software-card`、`.software-gallery`、`.software-tag` 等）
- 評論區段（`.press-source-title`、`.press-heading`、`.press-date`、`.press-author`）
- 語言切換規則（`.zh` / `.en` class）
- 回到頂部按鈕（`.back-to-top-shared`）
- 響應式斷點：手機 ≤768px、平板 769–1024px

#### `stats-data.js`
28 張學術數據卡片的中英文資料：

| 陣列 | 用途 | 更新時機 |
|------|------|---------|
| `statsZh[]` / `statsEn[]` | 學術數據卡片 | 論文/計畫數量變動時 |

> **維護提示**：日後更新卡片數據只需修改此檔，不需異動 `index.html` 或 `shared.js`。

#### `papers-data.js`
171 篇著作列表（5 大類別，含期別分組），中英文標題：

| 類別 | 篇數 |
|------|------|
| SCI / SSCI / EI 期刊論文 | 44 |
| 國際會議論文 | 49 |
| 專章、學位及其它論文 | 10 |
| 國內外媒體投書及採訪 | 49 |
| 國內研討會論文 | 19 |

> **維護提示**：新增或修改論文只需編輯此檔，然後執行 `python3 build-prerender.py` 重新產生預渲染 HTML。

#### `projects-data.js`
37 件研究計畫列表（3 種角色分類），含計畫名稱、期間、經費來源、編號、金額：

| 角色 | 件數 |
|------|------|
| 主持人 | 22 |
| 共同主持人 | 14 |
| 研究人員 | 1 |

> **維護提示**：新增或修改計畫只需編輯此檔，然後執行 `python3 build-prerender.py` 重新產生預渲染 HTML。

#### `build-prerender.py`
Python 建構工具，解析 `papers-data.js` 和 `projects-data.js`，將論文與計畫資料預渲染為靜態 HTML 並插入 `index.html`，同時產生 JSON-LD 結構化資料（ScholarlyArticle + ResearchProject），讓所有搜尋引擎（包括不執行 JavaScript 的爬蟲）都能索引完整內容。

```bash
# 更新論文或計畫資料後，執行以下指令重新產生預渲染 HTML
python3 build-prerender.py
```

#### `shared.js`
共用的 JavaScript 渲染函式與資料：

**資料陣列（修改數據請改這裡）：**
| 陣列 | 用途 | 更新時機 |
|------|------|---------|
| `honorsZh[]` / `honorsEn[]` | 榮譽列表 | 獲得新獎項時 |
| `alumniArticleZh[]` / `alumniArticleEn[]` | 老松國小傑出校友推薦文 | 推薦文內容變更時 |
| `honorsArticleZh[]` / `honorsArticleEn[]` | 杏壇芬芳錄推薦文 | 推薦文內容變更時 |
| `bioZh[]` / `bioEn[]` | 簡短自傳 | 自傳內容更新時 |
| `navItemsZh[]` / `navItemsEn[]` | 導覽列項目 | 新增/修改頁面連結時 |
| `footerLinks[]` | Footer logo 連結 | 新增/修改連結時 |

**核心函式：**
| 函式 | 用途 |
|------|------|
| `toggleLanguage()` | 中英文即時切換 |
| `showSection(id, event)` | 切換頁面內區段（介紹 / 榮譽 / 綜覽 / 著作 / 計畫 / 程式 / 評論） |
| `showEmbed(url, event)` | 以 iframe 嵌入顯示外部頁面（履歷 PDF） |
| `renderStats()` | 渲染學術數據卡片 |
| `renderHonors()` | 渲染榮譽列表（含照片/影片彈窗按鈕） |
| `showImageModal(src)` | 以彈窗顯示榮譽活動照片（`images/honors/` 目錄） |
| `showVideoModal(youtubeId)` | 以彈窗播放 YouTube 影片（如傑出校友專訪） |
| `renderAlumniArticle()` | 渲染老松國小傑出校友推薦文 |
| `renderHonorsArticle()` | 渲染杏壇芬芳錄推薦文（含 4 張照片） |
| `renderBio()` | 渲染簡短自傳 |
| `renderPapers()` | 渲染著作列表（讀取 `papers-data.js`） |
| `renderProjects()` | 渲染計畫列表（讀取 `projects-data.js`） |
| `renderNav()` | 渲染導覽列（高亮當前區段） |
| `renderFooter()` | 渲染 Footer logo 及提示框（使用縮圖） |
| `renderBackToTop()` | 初始化回到頂部按鈕（捲動超過 300px 時顯示） |
| `updateAriaLabels()` | 更新動態 aria-label（語言切換時同步更新） |
| `initShared()` | 初始化入口（無參數，自動偵測頁面） |

### 語言切換機制
- 使用 `.zh` / `.en` CSS class 控制元素顯示隱藏
- `html[lang="en"]` 時 `.en` 顯示、`.zh` 隱藏；預設中文時相反
- 切換時同步重新渲染 stats、honors、alumniArticle、honorsArticle、bio、papers、projects、nav（因資料陣列有中英文版本）

### HTML 區段順序

`index.html` 中 `<section>` 的排列順序與導覽列順序一致：

| 順序 | 區段 ID | 導覽列標籤（中/英） | 說明 |
|------|---------|-------------------|------|
| 1 | `bio` | 介紹 / About | 簡短自傳（預設顯示） |
| 2 | `overview` | 榮譽 / Honors | 責任與肯定 + 推薦文 |
| 3 | `lab` | 綜覽 / Overview | 28 張學術數據卡片 |
| 4 | `papers` | 著作 / Papers | 著作目錄 |
| 5 | `projects` | 計畫 / Projects | 計畫列表 |
| 6 | `software` | 程式 / Softwares | 開源程式 |
| 7 | `press` | 評論 / Op-Eds | 媒體投書（13 篇） |
| 8 | `embed` | 履歷 / CV | 內嵌 PDF 檢視 |

### 新增內容區段步驟
1. 在 `index.html` 的 `<main>` 中新增 `<section id="新區段" class="section">`（注意：須按導覽列順序插入正確位置）
2. 在 `shared.js` 的 `navItemsZh[]` 和 `navItemsEn[]` 加入新項目（設定 `section: '新區段'`）
3. 如需額外樣式，在 `shared.css` 中新增
4. 如有資料需預渲染，更新 `build-prerender.py` 並執行
5. 更新 `sitemap.xml` 的 `index.html` lastmod 日期

## 部署與發佈流程

本網站使用 GitHub Pages 自動部署，流程如下：

```
本機修改 → git add → git commit → git push origin master
                                        ↓
                              GitHub Pages 自動建置（約 30–60 秒）
                                        ↓
                              https://canslab1.github.io/ 生效
                                        ↓
                              IndexNow workflow 自動通知搜尋引擎
```

| 項目 | 說明 |
|------|------|
| 分支 | `master`（預設部署分支） |
| 建置方式 | GitHub Pages 內建（無需 Jekyll 或其他建置工具） |
| 生效時間 | push 後約 30–60 秒 |
| 快取 | GitHub CDN 快取約 10 分鐘，若未更新可嘗試 hard refresh（Ctrl+Shift+R） |
| 自訂域名 | 未設定（使用預設 `canslab1.github.io`） |
| HTTPS | GitHub Pages 自動提供 |

> **提示**：可在 GitHub → Settings → Pages 查看部署狀態，或透過 `gh run list` 查看 Actions 紀錄。

## 瀏覽器相容性

| 瀏覽器 | 最低版本 | 備註 |
|--------|---------|------|
| Chrome | 80+ | 完整支援 |
| Firefox | 78+ | 完整支援 |
| Safari | 14+ | 完整支援 |
| Edge | 80+ | 完整支援（Chromium 核心） |
| IE | ❌ 不支援 | 使用 ES6 語法及 CSS Grid |
| iOS Safari | 14+ | 響應式佈局已測試 |
| Android Chrome | 80+ | 響應式佈局已測試 |

### 響應式設計
- **手機**（≤768px）：單欄排列，卡片堆疊
- **平板**（769–1024px）：雙欄排列
- **桌面**（>1024px）：多欄排列，完整導覽列

### 效能考量
- 無框架依賴，純 HTML/CSS/JS，首次載入快速
- Footer logo 使用縮圖（`*-thumb.*`），減少約 1.1MB 頁面載入大小
- 圖片建議壓縮至 200KB 以下（詳見「圖片管理」章節）
- 外部資源：Ahrefs Analytics、Microsoft Clarity（非同步載入，不阻塞渲染）

## 家族故事頁面（stories.html）

| 項目 | 說明 |
|------|------|
| 路徑 | `stories.html` |
| 副標題 | 大頭仔的 60 道手工伍仁餡 |
| 用途 | 黃崇源教授親筆撰寫的 59 篇家族故事集 |
| 內容主題 | 從台北龍山寺旁的西園路到福建福州祖地，描繪三代人跨越時代的生活軌跡與情感連結 |
| 專屬樣式 | `stories.css`（獨立 CSS，不共用 `shared.css`） |
| 專屬腳本 | `stories.js`（獨立 JS，不共用 `shared.js`） |
| JSON-LD | `CollectionPage` + `CreativeWorkSeries`，標註作者、出版者、所屬網站 |
| og:type | `website` |
| og:image | `images/IMG-2.jpg`（教授大頭照） |
| 目錄功能 | 可收合的故事目錄（TOC），使用 `aria-expanded` / `aria-controls` 控制展開與收合 |
| 回到頂部 | 紫底白色向上箭頭圓形按鈕（與 index.html 的綠底按鈕風格區別） |

> **注意**：此頁面為完全獨立頁面，故事內容直接寫在 HTML 中（非透過 JS 陣列渲染），不依賴 `shared.css` / `shared.js`。

## 404 錯誤頁面

| 項目 | 說明 |
|------|------|
| 路徑 | `404.html` |
| 用途 | GitHub Pages 自訂 404 錯誤頁面 |
| 設計 | 自包含頁面，所有 CSS 以 `<style>` 內嵌（不引用外部樣式表） |
| robots | `noindex, follow`（不索引但允許追蹤連結） |
| 元素 | 教授大頭照、404 錯誤碼、提示訊息 |
| 導覽連結 | 「回首頁」（`/`）、「閱讀家族故事」（`stories.html`） |
| 背景圖片 | Header 使用 Unsplash 外部圖片（數位科技主題） |
| 分析追蹤 | Microsoft Clarity（ID: `rzlnthqbys`），以行內腳本載入 |
| 字型 | `system-ui` 系統字型堆疊（無外部字型依賴） |
| 響應式 | ≤768px 自動調整標題與錯誤碼字級 |

> **注意**：404.html 不引用 `shared.css` / `shared.js`，修改共用樣式不會影響此頁面。

### 子 repo 404 頁面

8 個研究軟體 repo 各自擁有獨立的 `404.html`，風格與主站一致（綠色漸層），提供三個導覽按鈕（README / GitHub Repo / Lab Home）。

## SEO 設定

### 各頁面 Meta 總覽

| 頁面 | og:type | og:image | JSON-LD 類型 | robots |
|------|---------|----------|-------------|--------|
| `index.html` | `profile` | `IMG-2.jpg`（大頭照） | `WebSite` + `ProfilePage` → `Person` + `ResearchOrganization` + 8× `SoftwareSourceCode` + `ItemList`（42 ScholarlyArticle + 22 ResearchProject） | `index, follow` |
| `stories.html` | `website` | `IMG-2.jpg`（大頭照） | `CollectionPage` + `CreativeWorkSeries` | `index, follow` |
| `404.html` | —（未設定） | —（未設定） | —（無） | `noindex, follow` |

### 共通 SEO 設定（所有可索引頁面）
- `og:locale`：`zh_TW`（主語言）
- `og:locale:alternate`：`en_US`（替代語言）
- `og:site_name`：`黃崇源教授個人網站`
- `og:video`：YouTube 傑出校友專訪嵌入（`index.html`）
- `twitter:card`：`summary_large_image`
- `hreflang`：`zh-TW`、`en`、`x-default` 三組互指
- `canonical`：各頁面指向自身完整 URL
- Schema.org：`<html>` 標籤加上 `itemscope itemtype="https://schema.org/WebPage"`
- `theme-color`：`#27ae60`（綠色主題色）
- `manifest.json`：PWA Web App Manifest（雙語描述）
- `apple-touch-icon`：iOS 桌面圖示（`images/icon-192.png`）
- `security.txt`：RFC 9116 安全聯絡資訊（根目錄）
- `feed.xml`：Atom feed（22 篇媒體投書 / 文章動態）

### 結構化資料（JSON-LD）詳細

**index.html**（整合八區段 + 內嵌 PDF）：
- `WebSite`：網站名稱、語言、描述
- `ProfilePage` → `Person`：姓名、別名（黃崇源、CY Huang、GSCOTT）、職稱、服務機構、學歷、研究領域、Google Scholar、ORCID
- `ResearchOrganization`：CANS Lab 組織資訊，隸屬長庚大學
- 8 個 `SoftwareSourceCode`：每套軟體的名稱、描述、程式語言、授權、GitHub 連結
- `ItemList`：42 篇 `ScholarlyArticle`（SCI/SSCI/EI 期刊論文）+ 22 個 `ResearchProject`（主持人計畫）— 由 `build-prerender.py` 自動產生

**stories.html**：
- `CollectionPage`：故事合集頁面
- `CreativeWorkSeries`：59 篇家族故事系列的名稱、描述、類型、作者

## LLM 搜尋設定

### llms.txt
根目錄下的 `llms.txt` 提供結構化純文字摘要，供 AI 爬蟲快速擷取教授資訊（姓名、職稱、研究領域、學術指標、獲獎、連結）。

### robots.txt
明確允許以下 AI 爬蟲存取：
- `GPTBot`、`ChatGPT-User`（OpenAI）
- `Google-Extended`（Google AI）
- `Claude-Web`、`Anthropic-AI`（Anthropic）
- `PerplexityBot`（Perplexity）
- `Bytespider`（ByteDance）
- `CCBot`（Common Crawl）
- `cohere-ai`（Cohere）
- `Applebot-Extended`（Apple Intelligence）
- `Meta-ExternalAgent`（Meta AI）
- `Amazonbot`（Amazon）
- `YouBot`（You.com）
- `Baiduspider`（百度）

### 網域驗證檔
| 檔案 | 用途 |
|------|------|
| `google9ba7e621fca59b66.html` | Google Search Console 驗證 |
| `BingSiteAuth.xml` | Bing Webmaster Tools 驗證 |
| `ahrefs_b1f1573b...` | Ahrefs SEO 工具驗證 |
| `openai-domain-verification=dv-...` | OpenAI 網域驗證 |
| `d22a81b36ccb45e085fe6679a822df52.txt` | IndexNow API key 驗證 |

## IndexNow 即時索引

### 簡介
IndexNow 是一個開放協議，可即時通知搜尋引擎（Bing、Yandex、Naver 等）網站內容已更新，加速索引收錄。

### 設定檔案
| 檔案 | 用途 |
|------|------|
| `d22a81b36ccb45e085fe6679a822df52.txt` | API key 驗證檔（放在根目錄供搜尋引擎驗證） |
| `.github/workflows/indexnow.yml` | GitHub Actions 自動提交 workflow |

### 自動提交（GitHub Actions）

workflow 會在以下時機自動執行，**不需要手動操作**：

| 觸發條件 | 說明 |
|----------|------|
| `push` 到 master | 當 `.html`、`.css`、`.js`、`.pdf`、`llms.txt`、`sitemap.xml`、`images/`、`articles/` 有變更時自動提交 |
| 每週一 09:00（台灣時間） | 定期排程，確保搜尋引擎持續收錄最新內容 |
| 手動觸發 | 在 GitHub → Actions → IndexNow Submit → Run workflow |

> **注意**：push 時只有內容檔案（html/css/js/pdf/images/articles）變更才會觸發，修改 README.md 等非內容檔案不會觸發。

### 子 repo 的 IndexNow

8 個研究軟體 repo（EpiRank、MV17、CASMIM、HETA、HATA、BCAT、SRAC-Agent、AED2）也各自設定了 IndexNow workflow，共用相同的 API key。觸發條件：push 到 main 分支修改 `.html`、`.css`、`.js`、`.md`、`images/` 時，或每週一 09:00 定期提交。

### 手動提交

如需立即手動觸發，可到 GitHub → Actions → IndexNow Submit → Run workflow。

## 子 Repo 架構與連結關係

### 概覽

本網站（`canslab1.github.io`）與 8 個研究軟體子 repo 形成一套統一的 GitHub Pages 站群架構。子 repo 利用 GitHub Pages 的多 repo 部署機制，自動掛載在主網域的子路徑下：

| 子 Repo | GitHub Pages URL | GitHub Repo |
|---------|-----------------|-------------|
| AED2 | https://canslab1.github.io/AED2/ | [canslab1/AED2](https://github.com/canslab1/AED2) |
| BCAT | https://canslab1.github.io/BCAT/ | [canslab1/BCAT](https://github.com/canslab1/BCAT) |
| CASMIM | https://canslab1.github.io/CASMIM/ | [canslab1/CASMIM](https://github.com/canslab1/CASMIM) |
| EpiRank | https://canslab1.github.io/EpiRank/ | [canslab1/EpiRank](https://github.com/canslab1/EpiRank) |
| HATA | https://canslab1.github.io/HATA/ | [canslab1/HATA](https://github.com/canslab1/HATA) |
| HETA | https://canslab1.github.io/HETA/ | [canslab1/HETA](https://github.com/canslab1/HETA) |
| MV17 | https://canslab1.github.io/MV17/ | [canslab1/MV17](https://github.com/canslab1/MV17) |
| SRAC-Agent | https://canslab1.github.io/SRAC-Agent/ | [canslab1/SRAC-Agent](https://github.com/canslab1/SRAC-Agent) |

### 雙向連結架構

主站與子 repo 之間存在雙向連結關係：

```
┌──────────────────────────────────────────────────────────┐
│  canslab1.github.io（主站 index.html）                    │
│  ├─ 「程式」區段 → 8 個軟體卡片                              │
│  │   每張卡片包含：                                         │
│  │   ├─ GitHub Repo 連結 → github.com/canslab1/{REPO}     │
│  │   └─ README 連結 → canslab1.github.io/{REPO}/          │
│  └─ Footer → github.com/canslab1                          │
└──────────────┬───────────────────────────────────────────┘
               │ 連結至子 repo 頁面
               ▼
┌──────────────────────────────────────────────────────────┐
│  canslab1.github.io/{REPO}/（子 repo index.html）         │
│  ├─ Header「CANS Lab」logo → canslab1.github.io/（主站）  │
│  ├─ Header「Lab Home」連結 → canslab1.github.io/（主站）   │
│  ├─ Header「GitHub Repo」→ github.com/canslab1/{REPO}     │
│  └─ Footer → canslab1.github.io/（主站）                   │
└──────────────────────────────────────────────────────────┘
```

- **主站 → 子 repo**：首頁「程式」區段（`<section id="software">`）的 8 張軟體卡片，每張包含 README 頁面連結（`https://canslab1.github.io/{REPO}/`）及 GitHub repo 連結。
- **子 repo → 主站**：每個子 repo 的 `index.html` header 包含「CANS Lab」logo 和「Lab Home」導覽連結，均指回 `https://canslab1.github.io/`；footer 也連回主站。

### 共用基礎設施

為了避免 8 個子 repo 各自維護相同的 CSS 與 JS，共用檔案統一託管在主站 repo（`canslab1.github.io`）中，由子 repo 透過絕對 URL 引用：

| 共用檔案 | 託管位置 | 用途 |
|---------|---------|------|
| `css/project-page.css` | 主站 repo | 子 repo 專案頁面的完整樣式（header、markdown-body、footer、loading spinner、響應式斷點） |
| `js/readme-loader.js` | 主站 repo | 自動偵測 repo 名稱、fetch README.md、以 marked.js 渲染 Markdown、以 highlight.js 語法高亮 |

子 repo 的 `index.html` 引用方式：
```html
<!-- CSS：替代原本 49 行的內嵌 <style> -->
<link rel="stylesheet" href="https://canslab1.github.io/css/project-page.css">

<!-- JS：替代原本 24 行的內嵌 <script> -->
<script src="https://canslab1.github.io/js/readme-loader.js"></script>
```

#### `css/project-page.css`
子 repo 專案頁面共用 CSS，包含：
- Reset 樣式與 body 排版（flex column、min-height 100vh）
- Header：綠色漸層背景（`#14532d` → `#1a7745`）、sticky 定位、品牌 logo 與導覽列
- Markdown body：白底圓角卡片、標題樣式（h1–h6 綠色邊框）、表格、程式碼區塊、引用
- Footer：白底分隔線、綠色連結
- Loading spinner：旋轉動畫、置中提示文字
- 響應式斷點：≤768px 自動調整 padding 與字級

#### `js/readme-loader.js`
README 渲染腳本，功能：
1. 從 URL 路徑自動偵測 repo 名稱（`/RepoName/` → `RepoName`）
2. Fetch 同目錄下的 `./README.md`
3. 使用 marked.js（GFM 模式）渲染為 HTML
4. 使用 highlight.js 進行程式碼區塊語法高亮
5. 錯誤處理：顯示錯誤訊息並附上 GitHub repo 連結

> **維護提示**：修改 `css/project-page.css` 或 `js/readme-loader.js` 會同時影響所有 8 個子 repo 的專案頁面，請謹慎更新。修改後需 push 主站 repo 並等待 GitHub Pages 部署（約 30–60 秒），子 repo 頁面即可載入更新後的檔案。

### 子 repo 模板結構

每個子 repo 的 `index.html` 遵循統一的模板結構（約 122 行），僅 repo 特定的 metadata 不同：

```html
<!DOCTYPE html>
<html lang="en">
<head>
    <!-- 1. Google Analytics（G-MXPNPC63XE，所有頁面共用） -->
    <!-- 2. Meta 標籤：title、description、og、twitter、keywords、robots -->
    <!-- 3. JSON-LD 結構化資料：SoftwareSourceCode + WebPage -->
    <!-- 4. Dublin Core metadata -->
    <!-- 5. Highwire Press（Google Scholar）metadata -->
    <!-- 6. 共用 CSS：project-page.css -->
    <!-- 7. highlight.js 主題 CSS -->
</head>
<body>
    <header>     <!-- CANS Lab logo + Lab Home / GitHub Repo 導覽 -->
    <main>       <!-- markdown-body 容器，載入 README.md -->
    <footer>     <!-- CANS Lab 版權聲明 + 主站連結 -->
    <!-- marked.js + highlight.js CDN -->
    <!-- 共用 JS：readme-loader.js -->
</body>
</html>
```

各子 repo 的差異僅在：
- `<title>` 和各 meta 標籤中的軟體名稱與描述
- JSON-LD 中的 `SoftwareSourceCode` 資訊（名稱、描述、程式語言、GitHub URL）
- Dublin Core 的 `DC.title`、`DC.subject`、`DC.description`
- Header 中的 GitHub Repo 連結指向各自的 repo

### 分支慣例

| Repo | 預設分支 | 部署分支 |
|------|---------|---------|
| canslab1.github.io（主站） | `master` | `master` |
| AED2、BCAT、CASMIM、EpiRank、HATA、HETA、MV17、SRAC-Agent | `main` | `main` |

> **注意**：主站使用 `master` 分支，8 個子 repo 使用 `main` 分支。Git 操作時請留意分支名稱。

## 分析工具與搜尋引擎管理

### 流量與行為分析

| 工具 | 用途 | 設定位置 |
|------|------|---------|
| **Ahrefs Analytics** | 流量分析、SEO 排名追蹤 | `index.html` `<head>` 載入 `analytics.ahrefs.com/analytics.js` |
| **Microsoft Clarity** | 用戶行為熱力圖、點擊分析、錄影回放 | `shared.js` 底部初始化，ID: `rzlnthqbys` |

### 搜尋引擎管理台

| 平台 | 網址 | 用途 |
|------|------|------|
| **Google Search Console** | [search.google.com/search-console](https://search.google.com/search-console) | 查看 Google 索引狀態、搜尋排名、點擊率、提交 sitemap |
| **Bing Webmaster Tools** | [bing.com/webmasters](https://www.bing.com/webmasters) | 查看 Bing/Yahoo 索引狀態、提交 URL、SEO 報告 |
| **Ahrefs Webmaster Tools** | [ahrefs.com/webmaster-tools](https://ahrefs.com/webmaster-tools) | 反向連結分析、SEO 健康度檢查 |

**建議定期檢查項目：**
- 每月查看 Google Search Console 的「涵蓋範圍」確認所有頁面已被索引
- 每月查看「成效」報告了解搜尋排名與點擊趨勢
- 如發現頁面未被索引，使用「網址審查」工具手動提交

## 第三方外部資源

本網站引用的外部服務與資源如下，全部以非同步方式載入，不阻塞頁面渲染：

| 資源 | 類型 | 載入位置 | 說明 |
|------|------|---------|------|
| **Ahrefs Analytics** | 流量分析 | 各頁面 `<head>`（`analytics.ahrefs.com/analytics.js`） | SEO 排名追蹤與流量分析 |
| **Microsoft Clarity** | 行為分析 | `shared.js` 底部 / `404.html` 行內腳本（`clarity.ms`，ID: `rzlnthqbys`） | 熱力圖、點擊分析、錄影回放 |
| **Unsplash 圖片** | 背景圖片 | `404.html` header（`images.unsplash.com`） | 數位科技主題背景圖（僅 404 頁面使用） |
| **ORCiD logo** | 外部圖片 | Footer（`orcid.org/assets/vectors/orcid.logo.icon.svg`） | ORCiD 官方 SVG logo |

### 字型策略
- 全站使用 `system-ui` 系統字型堆疊，**不依賴外部字型 CDN**（如 Google Fonts）
- 字型堆疊：`system-ui, -apple-system, 'PingFang TC', 'Microsoft JhengHei', 'Noto Sans TC', 'Helvetica Neue', Arial, sans-serif`
- 優點：零額外 HTTP 請求、更快的首次繪製、自動適配使用者作業系統

### DNS 預連接
各頁面 `<head>` 中設定 DNS 預取與預連接，減少第三方資源延遲：
- `dns-prefetch`：`scholar.google.com`、`www.clarity.ms`、`analytics.ahrefs.com`
- `preconnect`：`images.unsplash.com`

## 無障礙設計（Accessibility）

### 目前已實作的無障礙特性

| 特性 | 使用頁面 | 說明 |
|------|---------|------|
| `role="region"` | index / stories | 主要內容區段標記為語意區域 |
| `role="menubar"` / `role="menuitem"` | index | 導覽列使用 `<ul>` + `<li>` 結構，配合 ARIA 角色 |
| `aria-labelledby` | index / stories | 區段以標題 ID 作為可存取名稱 |
| `aria-label` | index | 綜覽區段（無標題）使用 `aria-label="綜覽"`；Footer 導覽使用 `aria-label="合作單位連結"` |
| `aria-expanded` | stories | 故事目錄展開/收合狀態通知螢幕閱讀器 |
| `aria-controls` | stories | 目錄按鈕與目錄內容區域的關聯 |
| `lang="zh-TW"` | 所有頁面 | 宣告頁面主語言，輔助螢幕閱讀器發音 |
| `alt` 屬性 | 所有頁面 | 所有 `<img>` 皆提供描述性替代文字 |
| 語意標記 | 所有頁面 | 使用 `<header>`、`<main>`、`<footer>`、`<section>`、`<nav>` 語意元素 |
| 鍵盤操作 | 所有頁面 | 語言切換按鈕、導覽連結均可透過鍵盤 Tab 鍵操作 |
| 響應式設計 | 所有頁面 | 手機 / 平板 / 桌面三種斷點自動適配 |

### 改善建議（未來可強化項目）
- 加入 skip navigation 連結（跳過導覽列直達主要內容）
- 為語言切換按鈕加入 `aria-label`（例如 `aria-label="切換至英文"`）
- 確保色彩對比度符合 WCAG 2.1 AA 標準（對比度 ≥4.5:1）

## 網站地圖（sitemap.xml）

| 頁面 | 優先度 | 更新頻率 | 備註 |
|------|--------|----------|------|
| `/` 和 `/index.html` | 1.0 | daily | 含 20 張 `images/honors/` 活動照片 |
| `/stories.html` | 0.8 | daily | |
| `/CV.pdf` | 0.7 | daily | |
| `/llms.txt` | 0.6 | daily | |
| `/feed.xml` | 0.5 | daily | Atom feed |
| `/security.txt` | 0.3 | yearly | RFC 9116 安全聯絡資訊 |
| `/manifest.json` | 0.3 | yearly | PWA Web App Manifest |
| `/humans.txt` | 0.3 | yearly | 網站製作者資訊 |
| `/CYHuang.html` | 0.1 | never | 重導向至 `index.html#papers` |
| `/software.html` | 0.1 | never | 重導向至 `index.html#software` |
| `/lab.html` | 0.1 | never | 重導向至 `index.html#lab` |

## 外部連結

| 平台 | 網址 |
|------|------|
| Google Scholar | https://scholar.google.com/citations?user=0klfzfAAAAAJ&hl=en |
| ORCID | https://orcid.org/0000-0002-8680-6755 |
| Chang Gung PURE | https://pure.lib.cgu.edu.tw/zh/persons/chung-yuan-huang-2 |
| Facebook | https://www.facebook.com/gscott.huang/ |
| Google Sites | https://sites.google.com/view/gscott-huang |
| GitHub | https://github.com/canslab1 |
| NSTC 學術研究 | https://arspb.nstc.gov.tw/（查詢黃崇源） |
| EpiRank（GitHub） | https://github.com/canslab1/EpiRank |
| EpiRank（README） | https://canslab1.github.io/EpiRank/ |
| MV17（GitHub） | https://github.com/canslab1/MV17 |
| MV17（README） | https://canslab1.github.io/MV17/ |
| CASMIM（GitHub） | https://github.com/canslab1/CASMIM |
| CASMIM（README） | https://canslab1.github.io/CASMIM/ |
| HETA（GitHub） | https://github.com/canslab1/HETA |
| HETA（README） | https://canslab1.github.io/HETA/ |
| HATA（GitHub） | https://github.com/canslab1/HATA |
| HATA（README） | https://canslab1.github.io/HATA/ |
| BCAT（GitHub） | https://github.com/canslab1/BCAT |
| BCAT（README） | https://canslab1.github.io/BCAT/ |
| SRAC-Agent（GitHub） | https://github.com/canslab1/SRAC-Agent |
| SRAC-Agent（README） | https://canslab1.github.io/SRAC-Agent/ |
| AED2（GitHub） | https://github.com/canslab1/AED2 |
| AED2（README） | https://canslab1.github.io/AED2/ |

## 內容管理指南

### 新增家族故事（stories.html）

`stories.html` 為獨立頁面，使用專屬的 `stories.css` 和 `stories.js`。

新增故事步驟：
1. 開啟 `stories.html`
2. 在故事列表區塊中，依照現有故事的 HTML 結構新增一筆
3. 故事內容直接寫在 HTML 中（非透過 JS 陣列渲染）
4. 如需新增圖片，放在 `images/` 目錄下
5. 更新 `sitemap.xml` 的 `stories.html` lastmod 日期

### 圖片管理

| 項目 | 建議 |
|------|------|
| 格式 | 照片優先使用 `.webp`（體積更小、品質佳）；logo / 圖示使用 `.png`（支援透明）；社群分享用 `og:image` 保留 `.jpg` 以確保平台相容性 |
| 尺寸 | 大頭照 / 活動照：最大寬度 800px；logo：最大寬度 200px；Footer 縮圖：寬度 128px 以內 |
| 檔案大小 | 單張建議 ≤200KB，大圖 ≤500KB |
| 命名規則 | 全小寫、用連字號分隔，例如 `honors1.png`、`epirank-screenshot.png`；縮圖加 `-thumb` 後綴 |
| 存放位置 | 一般圖片放 `images/`；榮譽活動照片放 `images/honors/`；軟體相關圖片放 `images/software/` |
| 壓縮工具 | 推薦 Python Pillow（`pip3 install Pillow`）轉 WebP；或 [TinyPNG](https://tinypng.com/)（線上）；macOS `sips` 不支援 WebP 寫入 |

**macOS 快速壓縮指令：**
```bash
# 將圖片寬度縮至 800px（等比縮放）
sips --resampleWidth 800 images/my-photo.jpg

# 製作 Footer 縮圖（寬度 128px）
sips -Z 128 images/logo.png --out images/logo-thumb.png
```

### 更新學術數據卡片

修改 `stats-data.js` 中的 `statsZh[]` 和 `statsEn[]` 陣列即可，每筆資料格式：
```javascript
{ number: '44', label: '期刊論文', url: 'https://...' }
```
- `number`：顯示的數字
- `label`：卡片標題
- `url`：點擊後前往的連結（可為空字串 `''`）

### 新增媒體投書（評論區段）

直接在 `index.html` 的 `<section id="press">` 中新增 `<li class="software-card">`，格式如下：
```html
<li class="software-card">
    <h4 class="press-heading">
        <span class="press-date">YYYY.MM.DD</span>
        <a href="文章網址" target="_blank" rel="noopener noreferrer">文章標題</a>
    </h4>
    <p class="press-author">作者名</p>
</li>
```
注意：同一報社內按日期由新到舊排列。

## 常見維護操作

| 操作 | 修改位置 |
|------|---------|
| 更新學術數據（論文數、引用數等） | `stats-data.js` → `statsZh[]` / `statsEn[]` |
| 新增/修改著作列表 | `papers-data.js` → 編輯後執行 `python3 build-prerender.py` |
| 新增/修改計畫列表 | `projects-data.js` → 編輯後執行 `python3 build-prerender.py` |
| 新增榮譽獎項 | `shared.js` → `honorsZh[]` / `honorsEn[]`（照片放 `images/honors/`，用 `showImageModal()` 按鈕；影片用 `showVideoModal()` 按鈕） |
| 更新傑出校友推薦文 | `shared.js` → `alumniArticleZh[]` / `alumniArticleEn[]` |
| 更新杏壇芬芳錄推薦文 | `shared.js` → `honorsArticleZh[]` / `honorsArticleEn[]` |
| 更新簡短自傳 | `shared.js` → `bioZh[]` / `bioEn[]` |
| 修改導覽列連結 | `shared.js` → `navItemsZh[]` / `navItemsEn[]` |
| 新增 Footer logo | `shared.js` → `footerLinks[]`（須同時製作縮圖） |
| 新增/修改媒體投書 | `index.html`「評論」區段（`<section id="press">`）的 `.software-card` 區塊 |
| 新增/修改研究軟體 | `index.html`「程式」區段（`<section id="software">`）中的 `.software-card` 區塊 |
| 更新履歷 | 替換 `CV.pdf` 檔案 |
| 修改頁面樣式 | `shared.css`（共用）或頁面內 `<style>`（專屬） |
| 更新 SEO 資訊 | `index.html` `<head>` 的 meta 標籤和 JSON-LD |
| 更新 LLM 摘要 | `llms.txt` |
| 更新網站地圖 | `sitemap.xml` |
| 通知搜尋引擎內容更新 | push 後由 GitHub Actions 自動執行 |
| 更換 Clarity 追蹤碼 | `shared.js` 底部 Clarity 初始化區塊 |

## 授權說明

本專案採用 [MIT License](LICENSE)，但以下資源**不包含**在 MIT 授權範圍內：

| 資源類型 | 授權狀態 |
|----------|---------|
| HTML / CSS / JavaScript 程式碼 | ✅ MIT 授權，可自由使用 |
| 研究軟體（EpiRank、HETA 等） | ✅ 各自 GitHub repo 採 MIT 授權 |
| 教授個人照片（`images/IMG-2.jpg`、`honors1–4.png`、`images/honors/*`） | ❌ 個人肖像權，未經授權不得使用 |
| `CV.pdf` 學術履歷 | ❌ 個人文件，僅供閱覽 |
| 各機構 logo（`csie.png`、`laosong-thumb.png` 等） | ❌ 各機構商標，僅限本網站使用 |
| 家族故事內容（`stories.html`） | ❌ 個人著作，版權所有 |

## 貢獻指南

本網站為個人學術網站，歡迎透過以下方式提供建議：

### 回報問題
- 前往 [GitHub Issues](https://github.com/canslab1/canslab1.github.io/issues) 建立新的 issue
- 請說明：問題描述、重現步驟、使用的瀏覽器與作業系統

### 提交修改
1. Fork 本專案
2. 建立新分支：`git checkout -b fix/issue-description`
3. 修改後提交：`git commit -m "Fix: 描述修改內容"`
4. 推送並建立 Pull Request

### 聯絡方式
| 管道 | 資訊 |
|------|------|
| 電子郵件 | gscott@mail.cgu.edu.tw |
| Facebook | [gscott.huang](https://www.facebook.com/gscott.huang/) |
| GitHub Issues | [canslab1.github.io/issues](https://github.com/canslab1/canslab1.github.io/issues) |
| 學校網頁 | [長庚大學資工系](https://csie.cgu.edu.tw/) |

## 備份與災難復原

### 備份策略

本網站以 Git 為核心版本控管，GitHub 遠端倉庫即為主要備份：

| 備份層級 | 說明 |
|---------|------|
| **Git 版本紀錄** | 每次 commit 都是完整快照，可隨時回溯任何歷史版本 |
| **GitHub 遠端倉庫** | `push` 後即自動同步至雲端，即使本機遺失仍可恢復 |
| **本機工作副本** | `/Users/gscott/Dropbox/...` 路徑下的工作目錄（Dropbox 自帶雲端同步） |

### 災難復原步驟

**情境一：本機檔案遺失**
```bash
# 重新從 GitHub 克隆
git clone https://github.com/canslab1/canslab1.github.io.git
```

**情境二：誤改檔案，想回到前一個版本**
```bash
# 查看某檔案的修改紀錄
git log --oneline -- 檔案名稱

# 恢復某檔案到特定 commit
git checkout <commit-hash> -- 檔案名稱

# 恢復後提交
git add 檔案名稱
git commit -m "Restore: 還原 檔案名稱 到 <commit-hash> 版本"
```

**情境三：想完整回退到某個 commit**
```bash
# 查看近期 commit
git log --oneline -10

# 建立回退 commit（安全方式，不改寫歷史）
git revert <commit-hash>
```

**情境四：GitHub Pages 部署異常**
1. 到 GitHub → Settings → Pages 確認部署分支為 `master`
2. 查看 Actions 面板是否有建置失敗
3. 如需強制重新部署，可做一個空 commit：
```bash
git commit --allow-empty -m "Trigger rebuild"
git push origin master
```

### 建議的額外備份措施
- 定期將 `CV.pdf` 等不在 Git 版控的重要檔案另行備份
- 考慮設定 GitHub 的 repository archive（Settings → Archive this repository）作為永久保存
- `images/` 目錄中的原始高解析圖片建議保留一份在本機或雲端硬碟

## 疑難排解（Troubleshooting）

### 語言切換無效
- **現象**：點擊中/英切換按鈕後頁面內容未改變
- **排查**：
  1. 開啟瀏覽器開發者工具（F12）→ Console 查看是否有 JS 錯誤
  2. 確認 `shared.js` 已正確載入（Network 面板）
  3. 確認 HTML 中 `.zh` / `.en` class 正確套用

### 學術數據卡片未顯示
- **現象**：點選「綜覽」區段後沒有顯示 28 張卡片
- **排查**：
  1. 確認 `stats-data.js` 在 `shared.js` 之前載入
  2. 確認 `stats-data.js` 語法正確（無多餘逗號或缺少引號）
  3. Console 是否出現 `ReferenceError: statsZh is not defined`

### 樣式破裂 / 排版異常
- **現象**：頁面排版混亂、卡片未對齊
- **排查**：
  1. 確認 `shared.css` 正確載入
  2. 嘗試 hard refresh：`Ctrl+Shift+R`（Windows）或 `Cmd+Shift+R`（Mac）
  3. 確認是否有 `<style>` 標籤衝突

### GitHub Pages 部署失敗
- **現象**：push 後網站未更新
- **排查**：
  1. 到 GitHub → Settings → Pages 確認部署分支為 `master`
  2. 查看 Actions 面板的 `pages build and deployment` 是否成功
  3. 等待 1–2 分鐘，GitHub Pages 需要時間建置
  4. 嘗試 hard refresh 清除瀏覽器快取

### IndexNow 提交失敗
- **現象**：GitHub Actions workflow 執行失敗
- **排查**：
  1. 到 GitHub → Actions → IndexNow Submit 查看失敗的 run log
  2. 確認 `d22a81b3...52.txt` 存在於網站根目錄
  3. 常見 HTTP 狀態碼：200/202（成功）、403（key 驗證失敗）、429（頻率限制，稍後再試）

### 論文/計畫預渲染失敗
- **現象**：執行 `python3 build-prerender.py` 後解析到 0 筆資料
- **排查**：
  1. 確認 `papers-data.js` / `projects-data.js` 語法正確（引號、逗號、括號配對）
  2. 確認 `index.html` 中有空的 `<div id="papers-container"></div>` 和 `<div id="projects-container"></div>`
  3. 如容器已含預渲染內容，需先還原為空再執行

### 外部連結重複開啟視窗
- **現象**：點擊「故事」「相簿」等導覽列外部連結時，每次都開新視窗
- **說明**：導覽列外部連結使用具名 `target`（如 `canslab-故事`），正常情況下重複點擊會跳到已開啟的視窗。但部分網站（如 Facebook、Google）設定了 `Cross-Origin-Opener-Policy: same-origin` 安全標頭，瀏覽器會清除具名視窗的關聯，導致每次都開新視窗。此為瀏覽器層級限制，無法透過 JavaScript 繞過。
