/* ===== Papers Data =====
 * 維護 index.html「著作」區段的論文列表
 * 修改本檔即可更新著作列表，不需異動 index.html 或 shared.js
 */

const papersData = [
    {
        titleZh: 'SCI / SSCI / EI 等級期刊論文',
        titleEn: 'SCI / SSCI / EI Journal Papers',
        noteZh: '註：Huang C.Y. 後面加上星號 * 表示本人是該篇期刊論文之通訊責任作者。',
        noteEn: 'Note: An asterisk (*) after Huang C.Y. indicates the corresponding author.',
        periods: [
            {
                titleZh: '長庚大學資訊工程學系正教授時期（2015.08 開始）',
                titleEn: 'Full Professor, Dept. of CSIE, Chang Gung University (since 2015.08)',
                items: [
            'Wang S.W., and Huang C.Y.* (2026) A Hybrid SVR-Based Framework for Cryptocurrency Price Forecasting and Strategy Backtesting. Applied Artificial Intelligence, 40(1), e2612793. https://doi.org/10.1080/08839514.2026.2612793. (SCI, IF: 4.3, Rank: 105/366 [27%] in Engineering, Electrical & Electronic)',
            'Wang S.W., Huang C.Y.*, and Sun C.T. (2022) Multiagent Diffusion and Opinion Dynamics Model Interaction Effects on Controversial Products. IEEE Access, 10, 115252-115270. https://doi.org/10.1109/ACCESS.2022.3218719. (SCI, IF: 4.098, Rank: 52/266 [20%] in Engineering, Electrical & Electronic)',
            'Huang C.Y. and Chin W.C.B. (2020) Distinguishing Arc Types to Understand Complex Network Strength Structures and Hierarchical Connectivity Patterns. IEEE Access, 8, 71021-71040. https://doi.org/10.1109/ACCESS.2020.2986017. (SCI, IF: 4.098, Rank: 52/266 [20%] in Engineering, Electrical & Electronic)',
            'Huang C.Y., Chin W.C.B., Fu Y.H., and Tsai Y.S. (2019) Beyond Bond Links in Complex Networks: Local Bridges, Global Bridges and Silk Links. Physica A: Statistical Mechanics and Its Applications, 536, 121027. https://doi.org/10.1016/j.physa.2019.04.263. (SCI, IF: 2.924, Rank: 26/85 [30%] in Physics, Multidisciplinary, EI)',
            'Huang C.Y., Chin W.C.B., Wen T.H., Fu Y.H., and Tsai Y.S. (2019) EpiRank: Modeling Bidirectional Disease Spread in Asymmetric Commuting Networks. Scientific Reports, 9:5415. https://doi.org/10.1038/s41598-019-41719-8. (SCI, IF: 4.122, Rank: 12/64 [19%] in Multidisciplinary Sciences)',
            'Fu Y.H., Huang C.Y.*, and Sun C.T. (2017) A Community Detection Algorithm Using Network Topologies and Rule-Based Hierarchical Arc-Merging Strategies. PLoS One, 12(11), e0187603. https://doi.org/10.1371/journal.pone.0187603. (SCI, IF: 2.806, Rank: 15/64 [24%] in Multidisciplinary Sciences)',
            'Lan Y.C., Wen T.H., Chang C.C., Liu H.F., Lee P.F., Huang C.Y., Chomel B.B., and Chen Y.M.A. (2017) Indigenous Wildlife Rabies in Taiwan: Ferret Badgers, A Long Term Terrestrial Reservoir. BioMed Research International, vol. 2017, Article ID 5491640, 6 pages. https://doi.org/10.1155/2017/5491640. (SCI, IF: 2.476 Rank: 67/160 [42%] in Biotechnology & Applied Microbiology)',
            'Fu Y.H., Huang C.Y.*, and Sun C.T. (2016) Using a Two-Phase Evolutionary Framework to Select Multiple Network Spreaders Based on Community Structure. Physica A: Statistical Mechanics and Its Applications, 461, 840-853. https://doi.org/10.1016/j.physa.2016.06.042. (SCI, IF: 2.243, Rank: 18/79 [23%] in Physics, Multidisciplinary, EI)',
            'Fu Y.H., Huang C.Y.*, and Sun C.T. (2015) Using Global Diversity and Local Topology Features to Identify Influential Network Spreaders. Physica A: Statistical Mechanics and Its Applications, 433, 344-355. https://doi.org/10.1016/j.physa.2015.03.042. (SCI, IF: 2.243, Rank: 18/79 [23%] in Physics, Multidisciplinary, EI)',
                ]
            },
            {
                titleZh: '長庚大學資訊工程學系副教授時期',
                titleEn: 'Associate Professor, Dept. of CSIE, Chang Gung University',
                items: [
            'Huang C.Y. (2015) An Agent-based Epidemic Simulation of Social Behaviors Affecting HIV Transmission among Taiwanese Homosexuals. Computational and Mathematical Methods in Medicine, vol. 2015, Article ID 867264, 10 pages. https://doi.org/10.1155/2015/867264. (SCI, IF: 0.887, Rank: 42/56 [75%] in Mathematical and Computational Biology)',
            'Lee C.L., Huang C.Y.*, Hsiao T.C., Wu C.Y., Chen Y.C, and Wang I.C. (2014) Impact of Vehicular Networks on Emergency Medical Services in Urban Areas. International Journal of Environmental Research and Public Health, 11(11), 11348-11370. https://doi.org/10.3390/ijerph111111348. (SCI, IF: 2.035, Rank: 101/225 [45%] in Environmental Sciences)',
            'Fu Y.H., Huang C.Y.*, and Sun C.T. (2014) Identifying Super-Spreader Nodes in Complex Networks. Mathematical Problems in Engineering, vol. 2014, Article ID 675713, 8 pages. https://doi.org/10.1155/2015/675713. (SCI, IF: 0.644, Rank 59/85 [70%] in Engineering, Multidisciplinary)',
            'Huang C.Y. and Wen T.H. (2014) Optimal Installation Locations for Automated External Defibrillators in Taipei 7-Elevens: Using GIS and a Genetic Algorithm with a New Stirring Operator. Computational and Mathematical Methods in Medicine, vol. 2014, Article ID 241435, 12 pages. https://doi.org/10.1155/2014/241435. (SCI, IF: 0.887, Rank: 42/56 [75%] in Mathematical and Computational Biology)',
            'Huang C.Y. and Lee C.L. (2014) Influences of Agents with a Self-Reputation Awareness Component in an Evolutionary Spatial IPD Game. PLoS One, 9(6), e99841. https://doi.org/10.1371/journal.pone.0099841. (SCI, IF: 3.057, Rank: 11/62 [18%] in Multidisciplinary Sciences)',
            'Huang C.Y. and Wen T.H. (2014) A Novel Private Attitude and Public Opinion Dynamics Model for Simulating Pluralistic Ignorance and Minority Influence. Journal of Artificial Societies and Social Simulation, 17(3), 8. https://doi.org/10.18564/jasss.2517. (SSCI, IF: 1.101, Rank: 28/93 [31%] in Social Sciences, Interdisciplinary)',
            'Wang S.W., Huang C.Y.*, and Sun C.T. (2014) Modeling Self-perception Agents in an Opinion Dynamics Propagation Society. Simulation: Transactions of the Society for Modeling and Simulation International, 90(3), 238-248. https://doi.org/10.1177/0037549713515029. (SCI, IF: 0.640, Rank: 77/106 [73%] in Computer Science, Software Engineering, EI)',
            'Wang S.W., Huang C.Y.*, and Sun C.T. (2014) A Proposal for Smart Agents with Adjusting Self-Reputation Capability for Preventing Fraud in Multi-Agent Societies. Journal of Information Science and Engineering, 30(2), 313-332. [PDF:papers/Journal%20Paper-17.pdf] (SCI, IF: 0.392, Rank: 133/143 [93%] in Computer Science, Information Systems)',
            'Huang C.Y.*, Wen T.H., and Tsai Y.S. (2013) FLUed: A Novel Four-Layer Model for Simulating Epidemic Dynamics and Assessing Intervention Policies. Journal of Applied Mathematics, vol. 2013, Article ID 325816, 20 pages. https://doi.org/10.1155/2013/325816. (SCI, IF: 0.834, Rank: 104/247 [43%] in Mathematics, Applied)',
            'Huang C.Y.*, Lee C.L., Wen T.H., and Sun C.T. (2013) A Computer Virus Spreading Model Based on Resource Limitations and Interaction Costs. The Journal of Systems and Software, 86(3), 801-808. https://doi.org/10.1016/j.jss.2012.11.027. (SCI, IF: 1.424, Rank: 31/105 [30%] in Computer Science, Theory and Methods)',
            'Huang C.Y.* and Sun C.T. (2012) Effects of Resource Limitations and Cost Influences on Computer Virus Epidemic Dynamics and Tipping Points. Discrete Dynamics in Nature and Society, vol. 2012, Article ID 473136, 15 pages. https://doi.org/10.1155/2012/473136. (SCI, IF: 0.646, Rank: 35/62 [57%] in Multidisciplinary Sciences)',
            'Tsai Y.S., Ko P.C.I., Huang C.Y., and Wen T.H. (2012) Optimizing Locations for the Installation of Automated External Defibrillators (AEDs) in Urban Public Streets Through the Use of Spatial and Temporal Weighting Schemes. Applied Geography, 35(1-2), 394-404. https://doi.org/10.1016/j.apgeog.2012.09.002. (SSCI, IF: 2.565, Rank: 10/77 [13%] in Geography)',
            'Huang C.Y.*, Tzou P.J., and Sun C.T. (2011) Collective Opinion and Attitude Dynamics Dependency on Informational and Normative Social Influences. Simulation: Transactions of the Society for Modeling and Simulation International, 87(10), 875-892. https://doi.org/10.1177/0037549710387940. (SCI, IF: 0.640, Rank: 77/106 [73%] in Computer Science, Software Engineering, EI)',
            'Huang C.Y.*, Wang S.W., and Sun C.T. (2011) Using Self-Aware Agents to Analyze Public Self-Consciousness in the Iterated Prisoner’s Dilemma. Simulation: Transactions of the Society for Modeling and Simulation International, 87(7), 600-615. https://doi.org/10.1177/0037549710391822. (SCI, IF: 0.640, Rank: 77/106 [73%] in Computer Science, Software Engineering, EI)',
            'Tsai Y.S., Huang C.Y., Wen T.H., Sun C.T., and Yen M.Y. (2011) Integrating Epidemic Dynamics with Daily Commuting Networks: Building a Multilayer Framework to Assess Influenza A (H1N1) Intervention Policies. Simulation: Transactions of the Society for Modeling and Simulation International, 87(5), 385-405. https://doi.org/10.1177/0037549710379481. (SCI, IF: 0.640, Rank: 77/106 [73%] in Computer Science, Software Engineering, EI)',
                ]
            },
            {
                titleZh: '長庚大學資訊工程學系助理教授兼資訊中心教學服務組組長時期',
                titleEn: 'Assistant Professor, Dept. of CSIE & Section Chief, IT Center, Chang Gung University',
                items: [
            'Tsai Y.S., Huang C.Y.*, and Sun C.T. (2011) Response to Wilson’s Note on “Influences of Resource Limitations and Transmission Costs on Epidemic Simulations and Critical Thresholds in Scale-Free Networks”. Simulation: Transactions of the Society for Modeling and Simulation International, 87(3), 267-270. https://doi.org/10.1177/0037549709351543. (SCI, IF: 0.640, Rank: 77/106 [73%] in Computer Science, Software Engineering, EI)',
            'Huang C.Y., Tsai Y.S., and Wen T.H. (2010) Simulations for Epidemiology and Public Health Education. Journal of Simulation. 4(1), 68-80. https://doi.org/10.1057/jos.2009.13. (SCI, IF: 1.164, Rank: 45/82 [55%] in Operations Research and Management Science, EI)',
            'Huang C.Y., Tsai Y.S., and Wen T.H. (2010) A Network-Based Simulation Architecture for Studying Epidemic Dynamics. Simulation: Transactions of the Society for Modeling and Simulation International. 86(5/6), 351-368. https://doi.org/10.1177/0037549709340733. (SCI, IF: 0.640, Rank: 77/106 [73%] in Computer Science, Software Engineering, EI)',
            'Tsai Y.S. and Huang C.Y.* (2010) Effects of Resources and Costs on Diffusion Dynamics. International Journal of Intelligent Information and Database Systems, 4(1), 19-42. https://doi.org/10.1504/IJIIDS.2010.031760. (EI)',
            'Huang C.Y. and Tsai Y.S. (2010) Effects of Friend-Making Resources/Costs and Remembering on Acquaintance Networks. Physica A: Statistical Mechanics and Its Applications, 389(3), 604-622. https://doi.org/10.1016/j.physa.2009.09.038. (SCI, IF: 1.772, Rank: 23/79 [30%] in Physics, Multidisciplinary, EI)',
            'Huang C.Y., Tsai Y.S., Sun C.T., Hsieh J.L., and Cheng C.Y. (2009) Influences of Resource Limitations and Transmission Costs on Epidemic Simulation and Critical Thresholds in Scale-Free Networks. Simulation: Transactions of the Society for Modeling and Simulation International, 85(3), 205-219. https://doi.org/10.1177/0037549709103775. (SCI, IF: 0.640, Rank: 77/106 [73%] in Computer Science, Software Engineering, EI)',
            'Hsieh J.L., Huang C.Y.*, Sun C.T., Tsai Y.S., and Kao G.Y.M. (2009) Learning to Build Network-Oriented Epidemic Simulation Models in Epidemiology Education. International Journal of Simulation and Process Modelling, 5(1), 31-41. https://doi.org/10.1504/ijspm.2009.025825. (EI)',
            'Huang C.Y. (2009) Teaching Epidemic Simulation and Policy-Making to Novice Researchers. International Journal of Instructional Media, 36(1), 107-117.',
            'Cheng C.Y., Huang C.Y.*, and Sun C.T. (2008) Mining Bridge and Brick Motifs from Complex Biological Networks for Functionally and Statistically Significant Discovery. IEEE Transactions on Systems, Man, and Cybernetics—Part B: Cybernetics, 38(1), 17-24. https://doi.org/10.1109/TSMCB.2007.908842. (SCI, IF: 2.361, Rank: 2/17 [12%] in Computer Science, Cybernetics, EI)',
            'Huang C.Y., Cheng C.Y., and Sun C.T. (2007) Bridge and Brick Network Motifs: Identifying Significant Building Blocks from Complex Biological Systems. Artificial Intelligence in Medicine, 41(2), 117-127. https://doi.org/10.1016/j.artmed.2007.07.006. (SCI, IF: 1.960, Rank: 28/94 [30%] in Computer Science, Artificial Intelligence, EI)',
            'Huang C.Y., Sun C.T., Cheng C.Y., and Hsieh J.L. (2007) Bridge and Brick Motifs in Complex Networks. Physica A: Statistical Mechanics and Its Applications, 377(1), 340-350. https://doi.org/10.1016/j.physa.2006.11.014. (SCI, IF: 1.772, Rank: 23/79 [30%] in Physics, Multidisciplinary, EI)',
            'Hsieh J.L., Sun C.T., Kao Y.M.G., and Huang C.Y.* (2006) Teaching through Simulation: Epidemic Dynamics and Public Health Policies. Simulation: Transactions of the Society for Modeling and Simulation International, 82(11), 731-759. https://doi.org/10.1177/0037549706074487. (SCI, IF: 0.640, Rank: 77/106 [73%] in Computer Science, Software Engineering, EI)',
            'Chen D.J., Tsai M.J., and Huang C.Y.* (2006) UI Design Patterns Generator for Pervasive Device. WSEAS Transactions on Computers, 5(9), 2114-2121. [PDF:papers/Journal%20Paper-37.pdf] (EI)',
            'Chen D.J., Tsai M.J., and Huang C.Y.* (2006) Visual Based Software Construction: Visual Requirement Authoring Tool and Visual Program Generator. WSEAS Transactions on Systems, 5(8), 1881-1888. [PDF:papers/Journal%20Paper-38.pdf] (EI)',
                ]
            },
            {
                titleZh: '元培科技大學資訊工程學系講師及助理教授兼系主任時期',
                titleEn: 'Lecturer & Assistant Professor / Dept. Chair, Dept. of CSIE, Yuanpei University of Medical Technology',
                items: [
            'Huang C.Y., Sun C.T., and Chu Y.W. (2006) A Co-adaptive Approach to XCS Classifier Systems. WSEAS Transactions on Systems, 5(2), 388-394. [PDF:papers/Journal%20Paper-39.pdf] (EI).',
            'Chu Y.W., Sun C.T., and Huang C.Y.* (2006) Regularity of Protein Secondary Structures and Its Prediction. WSEAS Transactions on Systems, 5(2), 380-385. [PDF:papers/Journal%20Paper-40.pdf] (EI).',
                ]
            },
            {
                titleZh: '交通大學資訊工程學系博士生暨中國科技大學資訊工程學系講師時期',
                titleEn: 'PhD Student, NCTU & Lecturer, China University of Technology',
                items: [
            'Huang C.Y., Hsieh J.L., Sun C.T., and Cheng C.Y. (2006) Teaching Epidemic and Public Health Policies through Simulation. WSEAS Transactions on Information Science and Applications, 3(5), 899-904. [PDF:papers/Journal%20Paper-41.pdf] (EI)',
            'Huang C.Y., Sun C.T., Hsieh J.L., Chen Y.M.A., and Lin H. (2005) A Novel Small-World Model: Using Social Mirror Identities for Epidemic Simulations. Simulation: Transactions of The Society for Modeling and Simulation International, 81(10), 671-699. https://doi.org/10.1177/0037549705061519. (SCI, IF: 0.640, Rank: 77/106 [73%] in Computer Science, Software Engineering, EI)',
            'Huang C.Y., Sun C.T., and Lin H.C. (2005) Influence of Local Information on Social Simulation in Small-World Network Models. Journal of Artificial Societies and Social Simulation, 8(4). https://www.jasss.org/8/4/8.html. (SSCI, IF: 1.101, Rank: 28/93 [31%] in Social Sciences, Interdisciplinary)',
            'Huang C.Y., Sun C.T., Hsieh J.L., and Lin H. (2004) Simulating SARS: Small-World Epidemiological Modeling and Public Health Policy Assessments. Journal of Artificial Societies and Social Simulation, 7(4). https://jasss.soc.surrey.ac.uk/7/4/2.html. (SSCI, IF: 1.101, Rank: 28/93 [31%] in Social Sciences, Interdisciplinary)',
                ]
            }
        ]
    },
    {
        titleZh: '國際會議論文',
        titleEn: 'International Conference Papers',
        periods: [
            {
                titleZh: null,
                titleEn: null,
                items: [
            'Huang C.Y. and Chin W.C.B. (2020) Hierarchical Arc Type Analysis (HATA) Algorithm: Incorporating Significant Direction Effects to Uncover Arc Strength in Complex Directed Networks. In Proceedings of the International School and Conference on Network Science (NetSci-X’20), Tokyo, Japan.',
            'Huang C.Y. and Chin W.C.B. (2018) Beyond Bond Links in Complex Networks: Local Bridges, Global Bridges and Silk Links. In Proceedings of the 6th International Conference on Information Technology: IoT and Smart City (ICIT’18), Hong Kong.',
            'Huang C.Y. (2018) Teaching Epidemic Policy-Making to Novice Researchers. In Proceedings of the 2018 International Conference on Education, Psychology, and Learning (ICEPL’18), Nagoya, Japan.',
            'Huang C.Y., Fu Y.H., and Sun C.T. (2018) A Fast Community Detection Algorithm. In Proceedings of the 20th International Conference on Computer Science, Cybersecurity and Information Technology (ICCSCIT’17), Tokyo, Japan.',
            'Huang C.Y. (2017) Influences of Self-Reputation Awareness Agents in an Evolutionary Spatial IPD Game. In Proceedings of the Tokyo 15th International Conference on Engineering & Technology, Computer, Basic & Applied Science (ECBA’17), Tokyo, Japan.',
            'Lee C.L., Huang C.Y., and Chen T.Y. (2017) Fast Conflict Detection for Two-Dimensional Packet Filters. In Proceedings of the 5th International Conference on Industrial Application Engineering (ICIAE’17), Kitakyushu, Japan.',
            'Fu Y.H., Huang C.Y.*, and Sun C.T. (2016) Using Network Topology and Rule-Based Strategy to Identify Community Structure in Social Networks. In Proceedings of the 4th International Conference on Industrial Application Engineering (ICIAE’16), Beppu, Japan.',
            'Huang C.Y. and Lee C.L. (2015) Self-Reputation Awareness Agents in Spatial IPD Game. In Proceedings of the 3rd International Conference on Industrial Application Engineering (ICIAE’15), Kitakyushu, Japan.',
            'Huang C.Y., Fu Y.H., and Sun C.T. (2014) Identify influential social network spreaders. In Proceedings of IEEE International Conference on Data Mining (ICDM’14), Shenzhen, China.',
            'Lee C.L., Huang C.Y., Hsiao T.C., Wu C.Y., Chen Y.C., and Wang I.C. (2014) Impact of Vehicular Ad Hoc Networks on Emergency Medical Services in Urban Areas. In Proceedings of 2014 Asia Global Land Project Conference (GLP’14), Taipei, Taiwan.',
            'Fu Y.H., Huang C.Y., and Sun C.T. (2014) Using Global Diversity and Local Features to Identify Influential Social Network Spreaders. In Proceedings of International Workshop on Curbing Collusive Cyber-gossips in Social Networks (C3’04) joint with IEEE/ACM International Conference on Advances in Social Networks Analysis and Mining (ASONAM’14), Beijing, China.',
            'Wen T.H. and Huang C.Y. (2014) Simulating Epidemic Dynamics and Assessing Intervention Policies. In Proceedings of International Conference on e-Commerce, e-Administration, e-Society, e-Education, and e-Technology (e-CASE & e-Tech’14), Nagoya, Japan.',
            'Wen T.H., Tsai Y.S., Ko P.C.I., and Huang C.Y. (2013) Using a Modified Genetic Algorithm in Identify Locations for Automated External Defibrillators. In Proceedings of International Conference on e-Commerce, e-Administration, e-Society, e-Education, and e-Technology (e-CASE & e-Tech’13), Kitakyushu, Japan. (Best Paper Award)',
            'Huang C.Y. (2012) Self-Aware Intelligent Agents in the Prisoner’s Dilemma. In Proceedings of Shanghai International Conference on Social Sciences (SICSS’12), Shanghai, China. (Session Chair)',
            'Wen T.H, Tsai Y.S., and Huang C.Y. (2012) Integration of Spatial and Temporal Variations of OHCA Patients and Emergency Service Facilities for Evaluating Potential Locations of Automated External Defibrillators (AED) in High Populated Urban Settings. In Proceedings of Global Change, Adaptation and Risk Management (GeoInformatics’12), Hong Kong.',
            'Huang C.Y., Tsai Y.S., and Wen T.H. (2011) Learning through Network-Based Simulations for Epidemiology Education. In Proceedings of Shanghai International Conference on Social Sciences (SICSS’11), Shanghai, China.',
            'Huang C.Y., Tsai Y.S., and Wen T.H. (2011) A Multilayer Framework to Assess Influenza Intervention Policies. In Proceedings of International Conference on Future Computer Sciences and Application (ICFCSA’11), Hong Kong, 132-136.',
            'Huang C.Y., Wang S.W., and Sun C.T. (2011) Self-Aware Intelligent Agents in the Prisoner’s Dilemma. In Proceedings of International Conference on Future Computer Sciences and Application (ICFCSA’11), Hong Kong, 127-131.',
            'Huang C.Y., Wang S.W., and Sun C.T. (2011) Modeling Agent Self-Awareness, Individual Performance and Collaborative Behavior. In Proceedings of IEEE 9th World Congress on Intelligent Control and Automation (WCICA’11), Taipei, Taiwan, 759-763.',
            'Huang C.Y. (2010) Influence of Self-Aware Agents with Public Self-Consciousness in the Iterated Prisoner’s Dilemma. In Proceedings of 9th WSEAS International Conference on System Science and Simulation in Engineering (ICOSSSE’10), Iwate, Japan, 133-138.',
            'Tsai Y.S., Huang C.Y., Wen T.H., and Sun C.T. (2010) A Multilayer Epidemiological Model Integrating Daily Commuting Network. In Proceedings of 9th WSEAS International Conference on System Science and Simulation in Engineering (ICOSSSE’10), Iwate, Japan, 77-83.',
            'Huang C.Y. and Wen T.H. (2010) A Multilayer Epidemic Simulation Framework Integrating Geographic Information System with Traveling Networks. In Proceedings of IEEE 8th World Congress on Intelligent Control and Automation (WCICA’10), Jinan, China, 2002-2007.',
            'Huang C.Y., Tsai Y.S., Wen T.H., and Sun C.T. (2010) A Multilayer Epidemic Simulation Framework Integrating Geographic Information System with Traveling Networks. In Proceedings of 29th IASTED International Conference on Modelling, Identification and Control (MIC’10), Innsbruck, Austria, 357-363.',
            'Wen T.Z., Tsai Y.S., Huang C.Y., and Yen M.Y. (2009) Multi-Layer Epidemic Dynamics Simulation (MEDSIM) for Spatial Transmission of Acute Respiratory Diseases. In Proceedings of 2009 Biosurveillance and Biosecurity Workshop International Perspectives, Taipei, Taiwan, R.O.C.',
            'Wang S.W., Sun C.T., and Huang C.Y. (2009) A Study of Agents with Self-awareness for Collaborative Behavior. In Proceedings of 14th Portuguese Conference on Artificial Intelligence (EPIA’09), Aveiro, Portugal, 425-435.',
            'Huang C.Y., Tsai Y.S., and Sun C.T. (2009) Effects of Resource and Remembering on Social Networks. In Proceedings of 8th International Conference on Autonomous Agents and Multiagent Systems (AAMAS’09), Budapest, Hungary, 837-842.',
            'Tsai Y.S., Sun C.T., and Huang C.Y. (2008) Epidemic Dynamics and Thresholds in Agent-Based Simulations under Realistic Resources and Cost Conditions. In Proceedings of IEEE/WIC/ACM International Conference on Intelligent Agent Technology (WI-IAT’08), Sydney, Australia, 65-70.',
            'Huang C.Y., Sun C.T., Cheng C.Y., and Tsai Y.S. (2008) Resource, Costs and Epidemic Critical Thresholds in Scale-Free Social Networks. In Proceedings of IEEE 7th World Congress on Intelligent Control and Automation (WCICA’08), Chongging, China, 1874-1879.',
            'Huang C.Y., Sun C.T., Cheng C.Y., and Tsai Y.S. (2008) Resource Limitations, Transmission Costs and Critical Thresholds in Scale-Free Networks. In Proceedings of 7th International Conference on Autonomous Agents and Multiagent Systems (AAMAS’08), Estoril, Portugal, 1121-1128.',
            'Huang C.Y., Hsieh J.L., Sun C.T., and Shih F.M. (2007) Mining Domain Ontological Information from Online Publications. In Proceedings of 6thd WSEAS International Conference on E-Activities (E-ACTIVITIES’07), Puerto De La Cruz, Tenerife, Canary Islands, Spain, 153-158.',
            'Cheng C.Y., Huang C.Y., and Sun C.T. (2007) Combining Motifs and Small-World Properties to Explore Complex Networks in Local and Global Views, In Proceedings of 2007 World Congress in Computer Science, Computer Engineering, and Applied Computing (WORLDCOMP’07), Las Vegas, Nevada, USA, 842-845.',
            'Huang C.Y., Cheng C.Y., and Sun C.T. (2007) Mining Bridge and Brick Network Motifs, In Proceedings of International Conference on Machine Learning and Data Mining (MLDM’07) Poster, Leipzig, Germany, 134-143.',
            'Sun C.T., Hsieh J.L., and Huang C.Y. (2006) Using Evolving Agents to Critique Subjective Music Compositions. In Proceedings of IEEE International Conference on Computational Intelligence and Security (CIS’06), Guangzhou, China, 474-480.',
            'Chen D.J., Tsai M.J., and Huang C.Y. (2006) Generating UI for Pervasive Devices Using Pattern-Based Approach. In Proceedings of 10th International Conference on Computers (CSCC’06), Athens, Greece, 837-843.',
            'Chen D.J., Tsai M.J., and Huang C.Y. (2006) Visual-Based Software Construction Methodology. In Proceedings of 10th International conference on Computers (CSCC’06), Athens, Greece, 831-836.',
            'Hsieh J.L., Sun C.T., and Huang C.Y. (2006) Implementing CASMIM for Epidemic Simulations. In Proceedings of First World Congress on Social Simulation (WCSS’06), Kyoto, Japan, 297-304.',
            'Huang C.Y., Sun C.T., Shih F.M., and Hsieh J.L. (2006) Mining Domain Ontological Information from Online Publications. In Proceedings of IADIS International Conference e-Society 2006, Dublin, Ireland, 433-440.',
            'Hsieh J.L., Sun C.T., and Huang C.Y. (2006) Using Evolving Agents to Critique Subjective Data: Recommending Music. In Proceedings of IEEE Congress on Evolutionary Computation (CEC’06), Vancouver, BC, Canada, 406-413.',
            'Huang C.Y., Hsieh J.L., Sun C.T., and Cheng C.Y. (2006) Teaching Epidemic Simulation and Policy-Making to Novice Researchers. In Proceedings of International Conference on Automatic Control, Modeling and Simulation (ACMOS’06), Prague, Czech Republic, 393-398.',
            'Cheng C.Y., Huang C.Y., Sun C.T., and Hsieh J.L. (2006) Bridge and Brick Network Motifs. In Proceedings of IEEE 6th World Congress on Intelligent Control and Automation (WCICA’06), Dalian, China, 1222-1226.',
            'Huang C.Y., Sun C.T., and Chu Y.W. (2005) Using a Coevolution Mechanism with a Dyna Architecture for Parameter Adaptation in XCS Classifier Systems. In Proceedings of International Conference on Computational Intelligence, Man-Machine Systems and Cybernetics (CIMMACS’05), Miami, Florida, USA, 173-178.',
            'Chu Y.W., Sun C.T., and Huang C.Y. (2005) Finding Regularity in Protein Secondary Structures using a Cluster-based Genetic Algorithm. In Proceedings of International Conference on Computational Intelligence, Man-Machine Systems and Cybernetics (CIMMACS’05), Miami, Florida, USA, 179-183.',
            'Sumodhee C., Hsieh J.L., Sun C.T., and Huang C.Y., and Chen Y.M.A. (2005) Impact of Social Behaviors on HIV Epidemic: A Computer Simulation View. In Proceedings of IEEE Computational Intelligence for Modelling Control and Automation (CIMCA’05), Vienna, Austria, 550-556.',
            'Hsieh J.L., Huang C.Y., Sun C.T., and Chen Y.M.A. (2005) Using the CAMIM Small-World Epidemic Model to Analyze Public Health Policies. In Proceedings of Western Multiconference on Health Sciences Simulation, New Orleans, LA, USA, 63-69.',
            'Huang C.Y. and Sun C.T. (2004) Self-Adaptive Routing Based on Learning Classifier Systems. In Proceedings of IEEE Congress on Evolutionary Computation (CEC’04), Portland, Oregon, USA, 678-682.',
            'Lin H.C., Huang C.Y., and Sun C.T. (2004) Influence of Local Information on Social Simulations under the Small-World Model, In Proceedings of IEEE Cyber Worlds, Tokyo, Japan, 356-360.',
            'Huang C.Y. and Sun C.T. (2004) Co-Adaptive Learning Classifier Systems Based on Coevolution within Dyna Architecture. In Proceedings of IEEE 5th World Congress on Intelligent Control and Automation (WCICA’04), Hangzhou, China, 2179-2183.',
            'Huang C.Y., Sun C.T., Yu C.H., and Chen C.F. (2004) Self-Adaptive Routing Method Based on Learning Classifier Systems. In Proceedings of IEEE 5th World Congress on Intelligent Control and Automation (WCICA’04), Hangzhou, China, 417-421.',
            'Yu C.H., Wu C.M., and Huang C.Y. (2004) A Study of Using Automatic Text Indexing to Analyze Web Browsing Behavior. In Proceedings of IEEE World Congress on Intelligent Control and Automation (WCICA’04), Hangzhou, China, 3991-3995.',
                ]
            }
        ]
    },
    {
        titleZh: '專章、學位及其它論文',
        titleEn: 'Book Chapters, Dissertations & Other Publications',
        periods: [
            {
                titleZh: null,
                titleEn: null,
                items: [
            'Wen T.H., Tsai Y.S., Huang C.Y., and Yen M.Y. (2009) Multi-Layer Epidemic Dynamics Simulation (MEDSim) for Spatial Transmission of Acute Respiratory Diseases. Intelligence and Security Informatics: Biosurveillance and Biosecurity, New York: West Publishing Co. (ISBN: 978-1-60317-007-9)',
            '黃崇源，孫春在，謝吉隆，林鶴玲（2009）第17章─模擬非典：基於小世界網絡的傳染病模型和公共衛生政策評估。聖昭瀚，張軍，杜建國（編輯）社會科學計算實驗理論與應用，中國上海：上海三聯出版社（南京大學工程管理學院研究所課程教科書，2012年江蘇省第十二屆哲學社會科學優秀成果一等獎，在環境行為與環境管理研究領域獲得國內同行的高度評價）。',
            'Huang C.Y., Cheng C.Y., and Sun C.T. (2008) Resource and Remembering Influences on Acquaintance Networks. Lecture Notes in Artificial Intelligence (LNAI), 4953, 281-291.',
            'Huang C.Y., Cheng C.Y., and Sun C.T. (2008) Resource Limitations, Transmission Costs and Critical Thresholds in Scale-Free Networks. Lecture Notes in Artificial Intelligence (LNAI), 4953, 485-494.',
            'Sun C.T., Hsieh J.L., and Huang C.Y.* (2007) Using Evolving agents to Critique Subjective Music Compositions. Lecture Notes in Artificial Intelligence (LNAI), 4456, 336-346.',
            'Huang C.Y., Hsieh J.L., and Sun C.T. (2006) Evaluating Subjective Compositions by the Cooperation between Human and Adaptive Agents. Lecture Notes in Artificial Intelligence (LNAI), 4293, 974-984. (SCIE, IF:0.513, Rank: 53/70 [76%] in Computer Science, Theory and Methods, EI)',
            'Huang C.Y. (2005) Small-World Epidemiological Modeling and Public Health Policy Assessment: Using the Social Mirror Identity Concept and Local Information for Network-based Epidemic Simulations: National Chiao Tung University, Hsinchu, Taiwan. (Ph.D. Dissertation)',
            'Huang C.Y. and Sun C.T. (2004) Parameter Adaptation within Co-adaptive Learning Classifier Systems. Lecture Notes in Computer Sciences (LNCS), 3103, 774-784.',
            '黃崇源，陳尚寬 (2002) 網際網路與資料庫的對話。青草湖社區大學E化工作坊，新竹市，台灣。',
            'Huang C.Y. (2000) Classifier System with Long Term Memory and Its Learning through Analogical Recognition: National Chiao Tung University, Hsinchu, Taiwan. (Master Thesis)',
                ]
            }
        ]
    },
    {
        titleZh: '國內外媒體投書及採訪',
        titleEn: 'Media Commentaries, Op-Eds & Press Coverage',
        periods: [
            {
                titleZh: null,
                titleEn: null,
                items: [
            'Huang C.Y. (2025) The beginning of your adventure with an AGI thinking machine! Taipei Times.',
            '黃崇源（2025）結交智能小幫手，結伴勇闖新世界。老松兒童45。',
            'Huang C.Y. (2025) Digging for gold or selling shovels? Taipei Times.',
            '黃崇源（2025）誰說賣鏟子的人知道金礦在那裡？自由時報，自由廣場。',
            '黃崇源（2024）從「老松國小校歌」談起。老松兒童43。',
            'Huang C.Y. (2023) Is ChatGPT actually intelligent? Taipei Times.',
            '黃崇源（2023）願景相承，生生不息。老松兒童42。',
            '黃崇源 (2021) 諾富特群聚事件升溫，學者揭兩大感染熱區。中國廣播公司新聞網。',
            '黃崇源 (2021) 提升疫苗接種率，學者籲政府擴大接種對象。中國廣播公司新聞網。',
            '黃崇源 (2021) 疫苗數量仍不足，學者憂無法堵缺口。中國廣播公司新聞網。',
            '黃崇源 (2021) 國內疫苗施打在即，學者呼籲確定兩件事。中國廣播公司新聞網。',
            '黃崇源 (2021) 學者研究疫苗施打率需7成以上，才有群體免疫。中國廣播公司新聞網。',
            '黃崇源 (2021) 避免衝擊雙北與竹科，學者籲桃園抗疫要守住。中國廣播公司新聞網。',
            'Huang C.Y. (2021) A message to students. Taipei Times.',
            '黃崇源 (2021) 流行病動態學教授給所有大學生的一封信。風傳媒。',
            'Huang C.Y. (2021) Taiwan’s pandemic history. Taipei Times.',
            '黃崇源（2021）一起來寫抗疫光榮史。自由時報，自由廣場。',
            'Huang C.Y. (2021) Getting the truth about virus shots. Taipei Times.',
            '黃崇源（2021）撥開副作用濃霧，正確看待接種計畫。自由時報，自由廣場。',
            '黃崇源 (2020) 機場「3高病毒量」設施！教授：候機室、貴賓室、公共廁所。東森新聞。',
            '黃崇源 (2020) 流行病學傳播預測，團隊：交通路徑關鍵。東森新聞。',
            '黃崇源 (2020) 資工學者疊合疾病分布，交通網正相關高達80%。TVBS新聞台。',
            '黃崇源 (2020) 通勤族注意！新冠肺炎傳播，專家列大台北「這5行政區」染疫高風險。中時新聞網。',
            '黃崇源 (2020) 台北市大安區、新北板橋、中和、三重和新莊，感染風險最高。藝術生活新聞網。',
            '黃崇源 (2020) 新冠病毒「通勤公司」！專家揭大台北5區高風險染疫。ETToday新聞雲。',
            '黃崇源 (2020) 通勤族當心！新冠病毒易散播，台北「5區」為染疫高風險。TVBS新聞台。',
            '黃崇源 (2020) 一屁股「坐走」病毒！學者：通勤、傳疫「高度」正相關。TVBS新聞台。',
            '黃崇源 (2020) 搭北捷、公車戴口罩，雙北五大紅區曝。三立新聞網。',
            '黃崇源 (2020) 病毒能「通勤」？專家曝大台北高風險區。中天新聞。',
            '黃崇源 (2020) 通勤易散播新冠病毒，專家曝台北5區風險高。中國時報。',
            '黃崇源 (2020) 雙北5大高風險區曝光，專家：通勤恐助長病毒傳播。今周刊。',
            '黃崇源 (2020) 病毒傳播恐與通勤高度相關！專家曝大台北5高風險區。自由時報。',
            '黃崇源 (2020) 國內新冠肺炎確診人數日日增，專家點出北北基5大高風險區。風傳媒。',
            '黃崇源 (2020) 通勤易傳播病毒？！雙北5區高風險區。華視新聞。',
            '黃崇源 (2020) 検疫指揮の「後藤新平はいないのか」台湾で揺らぐ「日本神話」新型肺炎対応に失望。日本產經新聞。',
            '黃崇源 (2020) 只囤泡麵衛生紙是錯的！學者籲政府教育正確儲量。聯合報。',
            '黃崇源 (2020) 學者研究：機場這3處，病毒量最高。聯合報。',
            '黃崇源 (2020) 學者研究：機場3處，病毒量最高。經濟日報。',
            '黃崇源 (2020) 機場3大地點，學者研究：病毒量最高。中時電子報。',
            '黃崇源 (2020) 機場「3大高病毒量設施」曝光！學者研究：連流感病毒也多。ETtoday新聞雲。',
            'Huang C.Y. (2020) Academics should help people, not sow fear. Taipei Times.',
            '黃崇源（2020）防疫，知識份子應有的責任絕非獵巫。風傳媒。',
            'Huang C.Y. and Chin W.C. (2020) Air travel significant to spread of COVID-19. Taipei Times.',
            '黃崇源，陳威全（2020）陳部長沒說的另一半答案。聯合報，民意論壇。',
            '黃崇源，陳威全（2020）塞翁失馬，焉知非福？台灣防疫做得最好！自由時報，自由共和國。',
            '黃崇源（2020）佛系防疫vs. 嚴謹防疫/日醫療優 台用閉關爭時間。聯合報，民意論壇。',
            'Huang C.Y. and Chin W.C. (2020) How to keep infection clusters from happening. Taipei Times.',
            '黃崇源，陳威全（2020）武漢肺炎進展到群聚、社區感染前，台灣可以怎麼做？自由時報，自由廣場。',
            '黃崇源 (2015) 李登輝只能代表自己。旺報。',
                ]
            }
        ]
    },
    {
        titleZh: '國內研討會論文',
        titleEn: 'Domestic Conference Papers',
        periods: [
            {
                titleZh: null,
                titleEn: null,
                items: [
            '陳禹杰，黃崇源（2021）使用流病風險指標評估新冠肺炎病毒封城策略成效，2021台灣網際網路研討會暨TANET & NCS全國計算機會議，東海大學，台中。',
            '黃崇源（2019）流行傳染病鄉鎮市流行評估指標，第17屆台塑企業應用技術研討會，長庚大學，桃園，台灣。',
            '黃崇源（2019）流行傳染病鄉鎮市流行評估指標，第17屆台塑企業應用技術研討會，長庚大學，桃園，台灣。',
            '黃崇源（2018）流行傳染病評估指標，第16屆台塑企業應用技術研討會，明志科技大學，台灣。',
            '賴誠信，黃崇源（2014）以基因演算法及時空間權重模式決定自動體外去顫器之台北市便利商店，2014年台灣地理資訊學會年會暨學術研討會，高雄，台灣。',
            '許仁賢，黃崇源（2014）基因演算法尋找SHA1碰撞的最佳解，2014年資訊技術應用及管理研討會，義守大學，高雄，台灣。',
            '賴誠信，黃崇源（2014）以基因演算法及時空間權重模式決定自動體外去顫器之台北市便利商店，2014年資訊技術應用及管理研討會，義守大學，高雄，台灣。',
            '黃崇源，陳思穎（2013）探討在社群網站上使用者隱私認知的不一致─以Facebook Applications為例，2013年資訊科技應用國際學術研討會，中國科技大學湖口校區，新竹，台灣。',
            '黃崇源，黃裕浩 (2012) 主動擴散式對等資源搜尋機制：透過資訊同步增進搜尋效率，第十七屆行動計算研討會，新北市三峽區大板根森林溫泉渡假村，新北市，台灣。',
            '黃崇源，黃逢傑（2012）以基因演算法結合人口統計資訊配置自動體外去顫器—以台北市為例，2012台灣地理資訊年會暨學術研討會，逢甲大學，台中，台灣。',
            '陳璽文，孫春在，黃崇源，溫在弘（2011）網路拓樸的遞移性：以傳染病的風險擴散為例，第二十八屆組和數學與計算理論研討會，澎湖科技大學，馬公市，澎湖，台灣。',
            '黃崇源，黃裕浩 (2008) 主動擴散式對等資源搜尋機制：透過資訊同步增進搜尋效率，TAAI 2008 第十三屆人工智慧與應用研討會，淡江大學蘭陽校區，宜蘭，台灣。',
            '鄭家胤，黃崇源，孫春在 (2008) 複雜網路之橋接式與磚式基調。生物資訊研討會，中華大學，新竹，台灣。',
            '黃崇源，鄭家胤 (2006) 橋接式與磚式網絡基調。2006中華決策科學學會年會暨學術研討會，新竹市，台灣，101頁。',
            '黃崇源，孫春在 (2006) 基於分類元系統之自適應路由路徑選擇機制。2006中華決策科學學會年會暨學術研討會，新竹市，台灣，102頁。',
            '黃崇源，孫春在 (2006) 共適應學習分類元系統。2006中華決策科學學會年會暨學術研討會，新竹市，台灣，103頁。',
            '黃崇源，孫春在，謝吉隆，陳宜民，林鶴玲 (2005) 小世界流行病學建模與公衛政策評估。健康與管理學術研討會論文摘要集，新竹市，台灣，87頁。',
            '孫春在，史馥銘，黃崇源 (2005) 從線上論文中擷取領域本體資訊。海峽二岸圖書館服務發展與創新高層論壇，北京，中國。',
            '黃崇源 (2001) 是否需要一份關於「電腦與網路」學程規劃指南？。全國社區大學教師研習會手冊，新竹市，台灣。',
                ]
            }
        ]
    }
];
