/* ===== Projects Data =====
 * 維護 index.html「計畫」區段的研究計畫列表
 * 修改本檔即可更新計畫列表，不需異動 index.html 或 shared.js
 */

const projectsData = [
    {
        titleZh: '主持人',
        titleEn: 'Principal Investigator',
        items: [
            { period: '109-110', source: '科技部專題研究計畫', grant: 'MOST-109-2314-B-182-037', internal: 'NERPD2K0481', amount: 'TWD$1,130,000', name: '以地理資訊系統為基礎之人工智慧基因演算法暨時空權重雙模型來探尋自動體外心臟除顫器於都會便利商店之最佳設置點' },
            { period: '107-108', source: '科技部專題研究計畫', grant: 'MOST-107-2221-E-182-069', internal: 'NERPD2H0551', amount: 'TWD$784,000', name: '整合資訊熵值及局域拓樸特性來識別複雜網絡超級傳播者' },
            { period: '106-107', source: '科技部專題研究計畫', grant: 'MOST-106-2221-E-182-073', internal: 'NERPD2G0431', amount: 'TWD$539,000', name: '使用規則式拓樸策略來快速識別複雜網絡社群及其階層結構' },
            { period: '105-106', source: '科技部專題研究計畫', grant: 'MOST-105-2221-E-182-069', internal: 'NERPD2F0341', amount: 'TWD$614,000', name: '線上熟人社群網絡演化之交友資源與記憶的影響' },
            { period: '104-105', source: '科技部專題研究計畫', grant: 'MOST-104-2221-E-182-047', internal: 'NERPD2E0351', amount: 'TWD$821,000', name: '反覆式囚犯困局中可自我聲譽調整能力之智慧型代理人的影響' },
            { period: '103-104', source: '科技部專題研究計畫', grant: 'MOST-103-2221-E-182-052', internal: 'NERPD2D0541', amount: 'TWD$500,000', name: '最適新型流感交通阻絕策略之基因演算法優選與防疫成本效益分析' },
            { period: '102-103', source: '國科會專題研究計畫', grant: 'NSC-102-2221-E-182-052', internal: 'NERPD2C0371', amount: 'TWD$761,000', name: '基於規範性與資訊性社會影響之多代理人意見與態度動態模擬' },
            { period: '101-102', source: '國科會專題研究計畫', grant: 'NSC-101-2119-M-182-001', internal: 'NERPD2B0301', amount: 'TWD$836,000', name: '整合型計畫之子計畫一：應用地理空間資料探礦與知識發掘技術於自願性地理資訊之擷取—以賞鳥資料為例' },
            { period: '100-101', source: '國科會專題研究計畫', grant: 'NSC-100-2221-E-182-058', internal: 'NERPD2A0471', amount: 'TWD$449,000', name: '公我意識之自我覺察代理人對於反覆囚犯困局之合作現象的影響與提升' },
            { period: '100-101', source: '新興病毒感染研究中心', grant: 'CMRPD290012', internal: null, amount: 'TWD$506,980', name: '台灣新型流感H1N1疫情之計量地理分析與時空傳播網絡監測平台 (2/2)' },
            { period: '99-100', source: '國科會專題研究計畫', grant: 'NSC-99-2314-B-182-031', internal: 'NERPD290581', amount: 'TWD$780,000', name: '建構「台灣通勤人口暨交通運輸網絡模型」探討新型流感傳播之最佳交通阻絕策略施行方式' },
            { period: '99-100', source: '新興病毒感染研究中心', grant: 'CMRPD290011', internal: null, amount: 'TWD$498,580', name: '台灣新型流感H1N1疫情之計量地理分析與時空傳播網絡監測平台 (1/2)' },
            { period: '98-99', source: '長庚大學校內研究計畫', grant: 'UERPD280351', internal: null, amount: 'TWD$497,500', name: '子計畫三：車用無線行動式智慧型P2P網路系統' },
            { period: '98-99', source: '新興病毒感染研究中心', grant: 'CMRPD260023', internal: null, amount: 'TWD$565,498', name: '利用社會網絡式傳染病模擬系統以協助新型流感監測並提供公衛參考 (3/3)' },
            { period: '98-99', source: '國科會專題研究計畫', grant: 'NSC-98-2314-B-182-043', internal: 'NERPD280521', amount: 'TWD$900,000', name: '台灣專屬新型流感傳播動態電腦模擬系統暨疫苗、抗病毒藥劑、居家隔離檢疫政策成本效益及最佳實施策略評估分析' },
            { period: '97-98', source: '國科會專題研究計畫', grant: 'NSC-97-2221-E-182-046', internal: 'NERPD270491', amount: 'TWD$570,000', name: '台灣重大流行性傳染病之公衛政策、疫苗、特效藥等防疫策略之流行病學電腦模擬與評估' },
            { period: '97-98', source: '長庚大學校內研究計畫', grant: 'UERPD270281', internal: null, amount: 'TWD$436,500', name: '主動擴散式P2P對等資源分享系統與搜尋演算法' },
            { period: '97-98', source: '新興病毒感染研究中心', grant: 'CMRPD260022', internal: null, amount: 'TWD$506,980', name: '利用社會網絡式傳染病模擬系統以協助新型流感監測並提供公衛參考 (2/3)' },
            { period: '96-97', source: '新興病毒感染研究中心', grant: 'CMRPD260021', internal: null, amount: 'TWD$578,580', name: '利用社會網絡式傳染病模擬系統以協助新型流感監測並提供公衛參考 (1/3)' },
            { period: '96-97', source: '國科會專題研究計畫', grant: 'NSC-96-2221-E-182-041', internal: 'NERPD260341', amount: 'TWD$558,000', name: '利用橋接式與磚式功能性基調分析探討實際複雜網絡的動態行為、設計原理與演化機制' },
            { period: '96-97', source: '疾病管制局研究計畫', grant: 'DOH96-DC-1102', internal: 'PERPD260011', amount: 'TWD$800,000', name: '臺灣地區腸病毒疫情之停課措施、洗手與戴口罩政策成效模擬與評估' },
            { period: '95-96', source: '國科會專題研究計畫', grant: 'NSC-95-2221-E-182-076', internal: 'NERPD250521', amount: 'TWD$449,000', name: '小世界流行病學建模與公共衛生政策之評估' }
        ]
    },
    {
        titleZh: '共同主持人',
        titleEn: 'Co-Principal Investigator',
        items: [
            { period: '109-110', source: '科技部專題研究計畫', grant: 'MOST-109-2221-E-992-069', internal: null, amount: null, name: '輿論擴散傳播模型運用人工智慧深度學習機制之探討' },
            { period: '108-111', source: '科技部專題研究計畫——沙克爾頓計畫（突破研究型）', grant: 'MOST-108-2368-H-002-MY2', internal: null, amount: null, name: '建立全球傳染疾病的傳播風險預警架構：評估東亞區域的境外移入風險' },
            { period: '102-103', source: '國科會專題研究計畫', grant: 'NSC-102-2410-H-002-168', internal: null, amount: null, name: '建立偵測時空群聚動態的分析方法：以追蹤傳染病擴散的傳播鏈為例' },
            { period: '102-103', source: '國科會專題研究計畫', grant: 'NSC-102-2221-E-182-034', internal: null, amount: null, name: '利用通用圖形處理器的多字串比對演算法之研究' },
            { period: '102-103', source: '國科會專題研究計畫', grant: 'NSC-102-2314-B-039-024', internal: null, amount: null, name: '台灣美沙酮病人自我管理即時互動手持行動系統' },
            { period: '102-102', source: '食品藥物管理局研究計畫', grant: '102TFDA-N-003', internal: null, amount: null, name: '結合網絡應用程式強化藥物濫用危害防制預警及宣導功能' },
            { period: '99-100', source: '國科會專題研究計畫', grant: 'NSC-99-2221-E-182-053', internal: null, amount: null, name: '點對點網路中檔案分享效率之研究' },
            { period: '98-99', source: '長庚大學校內研究計畫', grant: 'UERPD280331', internal: null, amount: null, name: '總計劃：行動式P2P無線感測網路之行車路徑動態導航系統' },
            { period: '98-99', source: '國科會專題研究計畫', grant: 'NSC-98-2221-E-182-033', internal: null, amount: null, name: '改善無結構式點對點網路服務品質之研究' },
            { period: '98-100', source: '國科會專題研究計畫', grant: 'NSC-98-2410-H-002-168-MY2', internal: null, amount: null, name: '台灣登革熱疫情的計量地理分析與空間模擬' },
            { period: '98-99', source: '疾病管制局研究計畫', grant: 'DOH98-NNB-1048', internal: null, amount: null, name: '利用風險分析方法建立藥物濫用預測模式及其防制教材' },
            { period: '98-99', source: '疾病管制局研究計畫', grant: 'DOH98-DC-1501', internal: null, amount: null, name: '台灣注射藥癮族群愛滋防治及減害計畫評估整合型計劃：探討影響藥癮愛滋減害計畫工作人員持續提供服務之相關因素調查' },
            { period: '97-98', source: '國科會專題研究計畫', grant: 'NSC-97-2314-B002-180', internal: null, amount: null, name: '建立登革熱疫情傳播網絡的時空風險模型' },
            { period: '95-96', source: '疾病管制局研究計畫', grant: 'DOH95-DC-1406', internal: null, amount: null, name: '利用整合式傳染病監測軟體以協助流行性感冒的監測' }
        ]
    },
    {
        titleZh: '研究人員',
        titleEn: 'Research Member',
        items: [
            { period: '98-99', source: '疾病管制局研究計畫', grant: 'DOH98-DC-1006', internal: null, amount: null, name: '疫災之風險管理及應變規劃研究' }
        ]
    }
];
