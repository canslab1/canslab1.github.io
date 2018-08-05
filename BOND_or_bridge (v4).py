# -*- encoding: utf-8 -*-

__author__   = '黃崇源, gscott@mail.cgu.edu.tw'
__date__     = '2012/03/01~2012/03/11'
__revision__ = '1.0'

import getopt
import math
import networkx as nx
import numpy
import os
import pylab
import sys
import shelve
import time
import scipy
import warnings
import xlwt
import scipy.cluster.hierarchy as hc

from matplotlib import pyplot as plot


"""
全域變數
"""
path                                = None                 # 目標網絡檔案的路徑
times                               = 100                  # 供目標網絡比較的隨機網絡的個數
debug                               = False                # 除錯模式，如果設定為 True，則會顯示訊息以告知使用者目前的執行狀況
quick                               = False
separation                          = 1

show_the_betweenness_result         = False
show_the_pagerank_result            = False
show_the_detailed_result            = True
show_the_major_result               = True
show_the_major_clustering_result    = True
show_the_detailed_clustering_result = True


"""
常數項，保持程式高可讀性且易於修改及維護
"""
STOP                                = 'stop'               # 停止下一階層的判斷
PASS                                = 'pass'               # 判斷尚未結束，下一階層繼續判斷

BOND                                = 'BOND'               # 鍵結連結，亦稱強連結，就算該連結消失，該連結兩端的節點還是能夠從共同的朋友，或者朋友的朋友，或者朋友的朋友的朋友
                                                           # (以此類推），獲得彼此的消息，或者傳遞消息
SINK                                = 'sink'               # 絲絮連結，一種橋接連結，但該連結某一端的節點的分支度為 1，表示此連結若消失，將造成該分支度為 1 的節點變成孤立節點
LOCAL_BRIDGE                        = 'local bridge'       # 區域橋接連結
GLOBAL_BRIDGE                       = 'global bridge'      # 全域橋接連結

SINK_COLOR                          = 'yellow'             # 目標網絡的絲絮連結的顏色
BOND_COLOR                          = 'blue'               # 目標網絡的鍵結連結的顏色
LOCAL_BRIDGE_COLOR                  = 'red'                # 目標網絡的區域橋接連結的顏色
GLOBAL_BRIDGE_COLOR                 = 'green'              # 目標網絡的全域橋接連結的顏色

SINK_BASIC_WIDTH                    = 1.0                  # 目標網絡的絲絮連結的基本寬度
BOND_BASIC_WIDTH                    = 1.0                  # 目標網絡的鍵結連結的基本寬度
BRIDGE_BASIC_WIDTH                  = 0.5                  # 目標網絡的橋接連結的基本寬度

NODE_SIZE_BASE                      = 140
NODE_SIZE                           =  80

REGULAR_NODE_COLOR                  = 'magenta'
IMPORTANT_NODE_COLOR                = 'pink'
SUPER_NODE_COLOR                    = 'red'

EGO_NETWORK                         = 'ego'


"""
演算法內部專用之常數項
"""
GRAPH_KEY_COMMON_NODES_LIST         = 'list_'
GRAPH_KEY_AVG_COMMON_NODES          = 'avg'
GRAPH_KEY_STD_COMMON_NODES          = 'std'
GRAPH_KEY_AVG_LIST                  = 'all.avg'
GRAPH_KEY_STD_LIST                  = 'all.std'
GRAPH_KEY_PASS_TO_NEXT_LAYER        = 'partial.w'
GRAPH_KEY_SHORTEST_PATH             = 'sp'
GRAPH_KEY_THRESHOLD_R1              = 'threshold.R1'
GRAPH_KEY_THRESHOLD_R2              = 'threshold.R2'

GRAPH_KEY_ENTROPY                   = 'entropy'
GRAPH_KEY_EDGE_CLASS                = 'edge_class'

GRAPH_KEY_NUMBER_OF_LAYER           = 'number_of_layer'

EDGE_KEY_LAYER                      = 'layer'
EDGE_KEY_COLOR                      = 'color'
EDGE_KEY_WIDTH                      = 'width'
EDGE_KEY_NEXT_STEP                  = 'next step'

NODE_KEY_EDGE_CLASS                 = 'edge_class'
NODE_KEY_NEW_ENTROPY                = 'new_entropy'
NODE_KEY_INFORMATION_GAIN           = 'information_gain'
NODE_KEY_GROUP_NUMBER               = 'group'


def debugmsg(s):
    """
    在除錯模式下顯示訊息以告知使用者目前的執行狀況
    
    參數：s 字串，顯示在 console 的字串
    """
    global debug
    if debug == True: print s


def generate_ego_graph(g, sp):
    for r in xrange(sp):
        for n in g.nodes_iter(data = False):
            if r == 0:
                g.node[n][EGO_NETWORK + str(r)] = set([n])
            else:
                g.node[n][EGO_NETWORK + str(r)] = set(g.node[n][EGO_NETWORK + str(r - 1)])
                for ng in nx.neighbors(g, n):
                    g.node[n][EGO_NETWORK + str(r)] = g.node[n][EGO_NETWORK + str(r)] | g.node[ng][EGO_NETWORK + str(r - 1)]


def get_ego_graph(g, s, t, l):
    index     = EGO_NETWORK + str(l - 1)
    node_list = set()
    for ng in nx.neighbors(g, s):
        if ng != t: node_list = node_list | g.node[ng][index]
    return node_list - set([s])


def compute_link_property(g, sp):
    """
    核心演算法：計算目標網絡的每一條連結的兩端節點在不同半徑下除了該連結之外的交集程度，以供稍後判斷 BOND/sink/local bridge/global bridge
    時間複雜度：O(m x l)
    m = 目標網絡的連結數目，在連通圖情況下，m 通常大於節點數目 n，卻遠小於節點數目的平方（n x n）
    l = 目標網絡的最短路徑，通常 l 遠小於 log(n)，可當作常數項 C 看待

    參數：g 目標網絡，必須是連通圖，若不是，函數 compute_link_property 將擷取目標網絡 g 最大的 component 來運算分析
    　　　sp 整數，通常是目標網絡的平均最短路徑，表示分析的階層數
    """
    c = g.copy()
    
    for i in xrange(sp): c.graph[GRAPH_KEY_COMMON_NODES_LIST + str(i + 1)] = []

    generate_ego_graph(c, sp)

    for s, t in g.edges_iter(data = False):
        #debugmsg('analyze the edge (' + s + ', ' + t + ')...')

        base_st_nodes = set([s, t])
        c.node[s][0]  = set()
        c.node[t][0]  = set()

        for i in xrange(sp):
            l = i + 1

            c.node[s][l] = get_ego_graph(c, s, t, l) - c.node[s][0] - base_st_nodes
            c.node[t][l] = get_ego_graph(c, t, s, l) - c.node[t][0] - base_st_nodes

            common_nodes = (c.node[s][l] & c.node[t][l]) | (c.node[s][l] & c.node[t][l - 1]) | (c.node[s][l - 1] & c.node[t][l])

            g[s][t][-l] = 0 if len(common_nodes) == 0 else float(len(common_nodes)) / (min(len(c.node[s][l]), len(c.node[t][l])) + min(len(c.node[s][l]), len(c.node[t][l - 1])) + min(len(c.node[s][l - 1]), len(c.node[t][l])))

            c.graph[GRAPH_KEY_COMMON_NODES_LIST + str(l)].append(g[s][t][-l])
                
            c.node[s][0] |= c.node[s][l]
            c.node[t][0] |= c.node[t][l]
                    
    for i in xrange(sp):
        l = str(i + 1)
        g.graph[GRAPH_KEY_AVG_COMMON_NODES + l] = scipy.mean(c.graph[GRAPH_KEY_COMMON_NODES_LIST + l])
        g.graph[GRAPH_KEY_STD_COMMON_NODES + l] = scipy.std( c.graph[GRAPH_KEY_COMMON_NODES_LIST + l])

    return g


def entropy (p):
    e = 0
    t = sum(p)
    for v in p:
        if not v == 0:
            pi = (float(v) / t)
            e += -(pi * math.log(pi, 2))
    return e


def network_clustering (g, layer):
    snapshot_g = {GLOBAL_BRIDGE:[], EDGE_KEY_LAYER + str(layer):[]}
    c          = g.copy()
    for s, t in g.edges_iter(data = False):
        if g[s][t][EDGE_KEY_LAYER + str(layer)].startswith(GLOBAL_BRIDGE):
            if c.has_edge(s, t): c.remove_edge(s, t)
        if g[s][t][EDGE_KEY_LAYER + str(layer)].startswith(SINK):
            if c.has_edge(s, t): c.remove_edge(s, t)
            tmpG = nx.Graph()
            if len(nx.neighbors(g, s)) == 0:
                if c.has_node(s): c.remove_node(s)
                g.node[s][NODE_KEY_GROUP_NUMBER] = "-0.01"
                tmpG.add_node(s)
            else:
                if c.has_node(t): c.remove_node(t)
                g.node[t][NODE_KEY_GROUP_NUMBER] = "-0.01"
                tmpG.add_node(t)
            snapshot_g[GLOBAL_BRIDGE].append(tmpG)
            snapshot_g[EDGE_KEY_LAYER + str(layer)].append(tmpG)
    snapshot_g[GLOBAL_BRIDGE].append(c)
    no = 1
    for sc in nx.connected_component_subgraphs(c):
        component_clustering(g, snapshot_g, sc, layer, "0." + ("%02d" % no))
        no += 1
    return snapshot_g
        

def component_clustering(bigG, sg, g, layer, cno):
    if g.order() == 1 or g.size() == 0 or layer == 0:
        for v in g.nodes_iter(data = False):
            bigG.node[v][NODE_KEY_GROUP_NUMBER] = cno
        if not layer == 0:
            sg[EDGE_KEY_LAYER + str(layer)].append(g)
        return
    c = g.copy()
    for s, t in g.edges_iter(data = False):
        if g[s][t][EDGE_KEY_LAYER + str(layer)] == LOCAL_BRIDGE + ' of layer ' + str(layer):
            if c.has_edge(s, t): c.remove_edge(s, t)
    sg[EDGE_KEY_LAYER + str(layer)].append(c)
    if layer > 1:
        sg[EDGE_KEY_LAYER + str(layer - 1)] = []
    no = 1
    for sc in nx.connected_component_subgraphs(c):
        component_clustering(bigG, sg, sc, layer - 1, cno + ("%02d" % no))
        no += 1                            


def link_analysis():
    """
    和隨機網絡相比，判斷目標網絡的每一條連結的拓樸連結屬性：BOND/sink/local bridge of layer #n/global bridge
    """
    global path, times, quick, separation

    # 剖析包含目錄及檔案名稱的 path 變數，並分別儲存成目錄 head、檔案主要名稱 tail 及副檔名 ext
    root, ext  = os.path.splitext(path)
    head, tail = os.path.split(root)

    # 假如目標網絡存在，則讀入記憶體，並計算其平均最短路徑，當作稍後計算基礎
    if os.path.exists(path) & os.path.isfile(path):
        debugmsg('read and analyse the target network...')
        # 打開 Pajek 網絡檔案，並轉換成無向圖
        G = nx.Graph(nx.read_pajek(path))
        # 設定開始的第一個 component graph 的編號為 0
        compNo = 0
        for g in nx.connected_component_subgraphs(G):
            # 如果這個 component 的節點數為 1 的話，則不做任何事
            if len(g.edges()) == 0: continue 

            # 計算平均最短路徑
            g.graph[GRAPH_KEY_SHORTEST_PATH] = nx.average_shortest_path_length(g)
            # component 的名稱等同 component 的編號
            g.name  = compNo
            compNo += 1

            # 決定每個節點要外看幾層，決定強弱連結
            layers = max(1, int(min((g.graph[GRAPH_KEY_SHORTEST_PATH] / 2.0), separation)) if quick else int(g.graph[GRAPH_KEY_SHORTEST_PATH] / 2.0))

            # 計算任意兩個有邊相連的節點的每一層共同朋友的正規化比率
            compute_link_property(g, layers)

            t_start = time.time()
            t_ttl   = 0
            rgs     = []
            # 產生供比較對應用的 times 個隨機網絡
            for c in xrange(times):
                rg_shelve_path = root + '_' + str(compNo) + '_' + str(c) + '_shelve.obj'
                rg_path        = root + '_' + str(compNo) + '_' + str(c) + ext
                rg_name        = 'random_network_' + str(g.name) + '_' + str(c)
                # 如果第 c 個隨機網絡過去已經產生過，則直接開檔讀取，否則重新建立一個隨機網絡
                if os.path.exists(rg_shelve_path) & os.path.isfile(rg_shelve_path):
                    debugmsg('read and analyse the random network #' + str(c) + ' from shelve file ' + rg_shelve_path + '...')                
                    sf = shelve.open(rg_shelve_path)
                    rg = sf[rg_name]
                else:
                    if os.path.exists(rg_path) & os.path.exists(rg_path):
                        debugmsg('read and analyse the random network #' + str(c) + ' from pajek file ' + rg_path + '...')
                        rg = compute_link_property(nx.connected_component_subgraphs(nx.Graph(nx.read_pajek(rg_path)))[0], layers)
                    else:
                        debugmsg('create, analyse and write the random network #' + str(c) + ' to pajek file ' + rg_path + '...')
                        rg = g.copy()
                        rg.graph['name'] = rg_name
                        if g.number_of_edges() > 2:
                            nx.connected_double_edge_swap(rg, g.number_of_edges())
                        compute_link_property(rg, layers)
                        nx.write_pajek(rg, rg_path)                    
                    rg.remove_nodes_from(rg.nodes())
                    sf = shelve.open(rg_shelve_path)
                    sf[rg_name] = rg
                sf.close()                    
                rgs.append(rg)
                t_inc  = time.time() - t_start
                t_ttl += t_inc
                debugmsg('+--- * Time spent (increment, total): (%f, %f)' % (t_inc, t_ttl))
                t_start = time.time()
            times = len(rgs)

            debugmsg('generate a threshold for BOND/bridge link analysis...')
            for i in xrange(layers):
                l = str(i + 1)
                g.graph[GRAPH_KEY_AVG_LIST + l] = []
                g.graph[GRAPH_KEY_STD_LIST + l] = []
                for j in xrange(times):
                    g.graph[GRAPH_KEY_AVG_LIST + l].append(rgs[j].graph[GRAPH_KEY_AVG_COMMON_NODES + l])
                    g.graph[GRAPH_KEY_STD_LIST + l].append(rgs[j].graph[GRAPH_KEY_STD_COMMON_NODES + l])
                g.graph[GRAPH_KEY_THRESHOLD_R1 + l] = scipy.mean(g.graph[GRAPH_KEY_AVG_LIST + l]) + 2 * scipy.mean(g.graph[GRAPH_KEY_STD_LIST + l])
                if g.graph[GRAPH_KEY_THRESHOLD_R1 + l] > 1:
                    g.graph[GRAPH_KEY_THRESHOLD_R1 + l] = 1.0

            debugmsg('assess the link property of every edge of the target network...')
            # phase 1: identify the sink links
            g.graph[SINK]          = 0
            g.graph[BOND]          = 0
            g.graph[LOCAL_BRIDGE]  = 0
            g.graph[GLOBAL_BRIDGE] = 0
            for s, t in g.edges_iter(data = False): 
                if (g.degree(s) == 1) | (g.degree(t) == 1):
                    g[s][t][EDGE_KEY_LAYER + '0'] = SINK
                    g[s][t][EDGE_KEY_NEXT_STEP]   = STOP
                    g[s][t][EDGE_KEY_WIDTH]       = SINK_BASIC_WIDTH
                    g[s][t][EDGE_KEY_COLOR]       = SINK_COLOR
                    g.graph[SINK] += 1
                else:
                    g[s][t][EDGE_KEY_NEXT_STEP]   = PASS

            # phase 2: identify the BOND/local bridge links on every layer
            for i in xrange(layers):
                l = -(i + 1)
                n = str(i + 1)            
                g.graph[GRAPH_KEY_PASS_TO_NEXT_LAYER + n] = []
                for s, t in g.edges_iter(data = False):                
                    if g[s][t][EDGE_KEY_NEXT_STEP] == STOP:
                        g[s][t][EDGE_KEY_LAYER + n] = g[s][t][EDGE_KEY_LAYER + str(i)]
                    elif g[s][t][l] >= g.graph[GRAPH_KEY_THRESHOLD_R1 + n]:                    
                        g[s][t][EDGE_KEY_LAYER + n] = BOND                        
                        g[s][t][EDGE_KEY_NEXT_STEP] = STOP
                        g[s][t][EDGE_KEY_WIDTH]     = (layers - i + 1) * BOND_BASIC_WIDTH
                        g[s][t][EDGE_KEY_COLOR]     = BOND_COLOR
                        g.graph[BOND] += 1
                    else:                        
                        g[s][t][EDGE_KEY_LAYER + n] = LOCAL_BRIDGE + ' of layer ' + n
                        g[s][t][EDGE_KEY_WIDTH]     = (layers - i + 1) * BRIDGE_BASIC_WIDTH
                        g[s][t][EDGE_KEY_COLOR]     = LOCAL_BRIDGE_COLOR
                        g.graph[GRAPH_KEY_PASS_TO_NEXT_LAYER + n].append(g[s][t][l])

                if len(g.graph[GRAPH_KEY_PASS_TO_NEXT_LAYER + n]) == 0:
                    g.graph[GRAPH_KEY_THRESHOLD_R2 + n] = 0
                else:
                    g.graph[GRAPH_KEY_THRESHOLD_R2 + n] = scipy.mean(g.graph[GRAPH_KEY_PASS_TO_NEXT_LAYER + n]) - scipy.std(g.graph[GRAPH_KEY_PASS_TO_NEXT_LAYER + n])
                    if g.graph[GRAPH_KEY_THRESHOLD_R2 + n] < 0:
                        g.graph[GRAPH_KEY_THRESHOLD_R2 + n] = 0.0
                    for s, t in g.edges_iter(data = False):
                        if g[s][t][EDGE_KEY_NEXT_STEP] == PASS:
                            if g[s][t][l] > g.graph[GRAPH_KEY_THRESHOLD_R2 + n]:                        
                                g[s][t][EDGE_KEY_NEXT_STEP] = STOP
                                g.graph[LOCAL_BRIDGE] += 1

            # phase 3: identify the global links
            for s, t in g.edges_iter(data = False):
                if g[s][t][EDGE_KEY_NEXT_STEP] == PASS:
                    g[s][t][EDGE_KEY_LAYER + n] = GLOBAL_BRIDGE
                    g[s][t][EDGE_KEY_WIDTH]     = BRIDGE_BASIC_WIDTH
                    g[s][t][EDGE_KEY_COLOR]     = GLOBAL_BRIDGE_COLOR
                    g.graph[GLOBAL_BRIDGE] += 1

            # extra phase 4: identify the node entropy
            ns = []
            nc = []
            g.graph[GRAPH_KEY_EDGE_CLASS] = {BOND:g.graph[BOND], LOCAL_BRIDGE:g.graph[LOCAL_BRIDGE], GLOBAL_BRIDGE:g.graph[GLOBAL_BRIDGE]}
            g.graph[GRAPH_KEY_ENTROPY]    = entropy(g.graph[GRAPH_KEY_EDGE_CLASS].values())
            for s in g.nodes_iter(data = False):
                g.node[s][NODE_KEY_EDGE_CLASS] = g.graph[GRAPH_KEY_EDGE_CLASS].copy()
                for t in nx.neighbors(g, s):
                    for key in g.node[s][NODE_KEY_EDGE_CLASS].keys():
                        if g[s][t][EDGE_KEY_LAYER + str(layers)].startswith(key):
                            g.node[s][NODE_KEY_EDGE_CLASS][key] -= 1
                g.node[s][NODE_KEY_NEW_ENTROPY]      = entropy(g.node[s][NODE_KEY_EDGE_CLASS].values())
                g.node[s][NODE_KEY_INFORMATION_GAIN] = max(0, g.graph[GRAPH_KEY_ENTROPY] - g.node[s][NODE_KEY_NEW_ENTROPY])
                ns.append(g.node[s][NODE_KEY_INFORMATION_GAIN])
                nc.append([REGULAR_NODE_COLOR, IMPORTANT_NODE_COLOR, SUPER_NODE_COLOR][max(0, int(math.ceil(g.node[s][NODE_KEY_INFORMATION_GAIN])))])
            ns_avg = scipy.mean(ns)
            if not ns_avg == 0:
                ns = [NODE_SIZE_BASE + NODE_SIZE * (value / ns_avg) for value in ns]

            # extra phase 5: save the network fingerprint into a file
            nfp_shelve_path = 'network_fingerprints.obj'
            if os.path.exists(nfp_shelve_path) & os.path.isfile(nfp_shelve_path):
                sf = shelve.open(nfp_shelve_path)
                finger_prints = sf['finger_prints']
            else:
                sf = shelve.open(nfp_shelve_path)                
                finger_prints = {}
            d = float(g.graph[BOND] + g.graph[LOCAL_BRIDGE] + g.graph[GLOBAL_BRIDGE] + g.graph[SINK])
            network_name = root + '_' + str(compNo)
            finger_prints[network_name] = {0:round(g.graph[BOND]          / d, 4),
                                           1:round(g.graph[LOCAL_BRIDGE]  / d, 4),
                                           2:round(g.graph[GLOBAL_BRIDGE] / d, 4),
                                           3:round(g.graph[SINK]          / d, 4)}
            corr_table = {}
            for net_name1, net_series1 in finger_prints.items():
                corr_table[net_name1] = {}
                for net_name2, net_series2 in finger_prints.items():
                    corr_table[net_name1][net_name2] = numpy.corrcoef(net_series1.values(), net_series2.values())[0, 1]
            sf['corr_table']    = corr_table
            sf['finger_prints'] = finger_prints
            sf.close()

            debugmsg('write the analysis results to a pajek file...')
            ng = nx.Graph()
            ng.add_nodes_from(g.nodes())
            ng.add_edges_from(g.edges(data = True))
            ng.graph['name'] = root + '_' + str(compNo) + '_result' + ext
            nx.write_pajek(ng, root + '_' + str(compNo) + '_result' + ext)

            debugmsg('write the analysis results to a excel file...')
            # Phase 1: write texphe analysis results of the target network to the sheet1
            c = g.copy()
            book = xlwt.Workbook()
            s1 = book.add_sheet('target network')
            s2 = book.add_sheet(str(times) + ' random networks')
            s3 = book.add_sheet('node information')
            si = xlwt.Style.easyxf('font: name Arial, colour dark_red, bold True; alignment: horizontal left;')
            st = xlwt.Style.easyxf('font: name Arial, colour dark_red, bold True; alignment: horizontal center;')
            sb = xlwt.Style.easyxf('font: name Arial, colour dark_blue;')

            # phase 1.1: create the heading data of the analysis report
            numpy.seterr(all='ignore')
            row = 5
            col = 7
            s1.write( 0, 0, 'number of nodes = ' + str(g.number_of_nodes()),                             si)
            s1.write( 1, 0, 'number of edges = ' + str(g.number_of_edges()),                             si)
            s1.write( 2, 0, 'average degree = '  + str(g.number_of_edges() * 2.0 / g.number_of_nodes()), si)
            s1.write( 3, 0, 'diameter = '        + str(nx.diameter(g)),                                  si)
            s1.write( 4, 0, 'average shortest path = '            + str(round(g.graph[GRAPH_KEY_SHORTEST_PATH], 4)),       si)
            s1.write( 5, 0, 'average clustering coefficient = '   + str(round(nx.average_clustering(g), 4)),               si)
            s1.write( 6, 0, 'degree assortativity coefficient = ' + str(round(nx.degree_assortativity_coefficient(g), 4)), si)
            s1.write( 7, 0, 'BOND = '          + str(g.graph[BOND])              + " (" + str(100 * round(float(g.graph[BOND])          / g.size(), 4)) + '%)', si)
            s1.write( 8, 0, 'sink = '          + str(g.graph[SINK])              + " (" + str(100 * round(float(g.graph[SINK])          / g.size(), 4)) + '%)', si)
            s1.write( 9, 0, 'local bridge = '  + str(g.graph[LOCAL_BRIDGE])      + " (" + str(100 * round(float(g.graph[LOCAL_BRIDGE])  / g.size(), 4)) + '%)', si)
            s1.write(10, 0, 'global bridge = ' + str(g.graph[GLOBAL_BRIDGE])     + " (" + str(100 * round(float(g.graph[GLOBAL_BRIDGE]) / g.size(), 4)) + '%)', si)
            s1.write(11, 0, 'entropy = '       + str(g.graph[GRAPH_KEY_ENTROPY]), si)
            s1.write(row - 1, col - 6, 'st.sp',  st)
            s1.write(row - 1, col - 5, 'avg.sp', st)
            s1.write(row - 1, col - 4, 's.cc',   st)
            s1.write(row - 1, col - 3, 't.cc',   st)
            s1.write(row - 1, col - 2, 'source', st)
            s1.write(row - 1, col - 1, 'target', st)        
            for i in xrange(layers):
                s1.write(row - 3, col + (i * 2),     'R1 = ' + str(round(g.graph[GRAPH_KEY_THRESHOLD_R1 + str(i + 1)], 4)), si)
                s1.write(row - 2, col + (i * 2),     'R2 = ' + str(round(g.graph[GRAPH_KEY_THRESHOLD_R2 + str(i + 1)], 4)), si)
                s1.write(row - 1, col + (i * 2),     'intersection weight', st)
                s1.write(row - 1, col + (i * 2) + 1, 'layer ' + str(i + 1), st)

            # phase 1.2: create the body data of the analysis report
            for s, t in g.edges_iter(data = False):
                """
                速度很慢，修改中。
                c.remove_edge(s, t)
                try:
                    st_sp     = nx.shortest_path_length(c, s, t)
                    st_avg_sp = round(nx.average_shortest_path_length(c), 2)
                except nx.NetworkXNoPath:
                    st_sp     = 'unlimited'
                    st_avg_sp = '2 component'

                s1.write(row, col - 6, st_sp,                         sb)
                s1.write(row, col - 5, st_avg_sp,                     sb)
                """
                s1.write(row, col - 6, "unused",                      sb)
                s1.write(row, col - 5, "unused",                      sb)
                s1.write(row, col - 4, round(nx.clustering(g, s), 2), sb)
                s1.write(row, col - 3, round(nx.clustering(g, t), 2), sb)
                s1.write(row, col - 2, s, sb)
                s1.write(row, col - 1, t, sb)
                """
                c.add_edge(s, t)
                """
                for i in xrange(layers):
                    s1.write(row, col + (i * 2), round(g[s][t][-(i + 1)], 4), sb)
                    if (i == 0):                    
                        s1.write(row, col + (i * 2) + 1, g[s][t][EDGE_KEY_LAYER + str(i + 1)], sb)
                    elif (g[s][t][EDGE_KEY_LAYER + str(i + 1)] != g[s][t][EDGE_KEY_LAYER + str(i)]):
                        s1.write(row, col + (i * 2) + 1, g[s][t][EDGE_KEY_LAYER + str(i + 1)], sb)
                    else:
                        s1.write(row, col + (i * 2) + 1, '...', sb)
                row += 1

            # phase 2: write the analysis results of the random networks to the sheet2
            row = 5
            col = 3
            for i in xrange(layers):
                l = str(i + 1)
                s2.write(row - 2, col + (i * 4),     'layer ' + l, st)
                s2.write(row - 1, col + (i * 4),     'AVG',        st)
                s2.write(row - 1, col + (i * 4) + 1, 'STD',        st)                     
                for j in xrange(times):
                    s2.write(row + j, col + (i * 4),     rgs[j].graph[GRAPH_KEY_AVG_COMMON_NODES + l], sb)
                    s2.write(row + j, col + (i * 4) + 1, rgs[j].graph[GRAPH_KEY_STD_COMMON_NODES + l], sb)

            # extra phase 3: write the analysis results of node entropy
            row = 1
            col = 1
            now = 1
            s3.write(row,     col + 0, 'node',                                       st)
            s3.write(row - 1, col + 1, 'o.entropy = ',                               st)            
            s3.write(row,     col + 1, 'n.entropy',                                  st)
            s3.write(row - 1, col + 2, g.graph[GRAPH_KEY_ENTROPY],                   sb)
            s3.write(row,     col + 2, 'gain',                                       st)
            s3.write(row - 1, col + 3, g.graph[GRAPH_KEY_EDGE_CLASS][BOND],          sb)
            s3.write(row,     col + 3, 'BOND',                                       st)
            s3.write(row - 1, col + 4, g.graph[GRAPH_KEY_EDGE_CLASS][LOCAL_BRIDGE],  sb)
            s3.write(row,     col + 4, 'local bridge',                               st)
            s3.write(row - 1, col + 5, g.graph[GRAPH_KEY_EDGE_CLASS][GLOBAL_BRIDGE], sb)
            s3.write(row,     col + 5, 'global bridge',                              st)
            s3.write(row,     col + 6, 'avg',                                        st)
            s3.write(row + 1, col + 6, ns_avg,                                       sb)
            for s in g.nodes_iter(data = False):
                s3.write(row + now, col + 0, s,                                             sb)
                s3.write(row + now, col + 1, g.node[s][NODE_KEY_NEW_ENTROPY],               sb)
                s3.write(row + now, col + 2, g.node[s][NODE_KEY_INFORMATION_GAIN],          sb)
                s3.write(row + now, col + 3, g.node[s][NODE_KEY_EDGE_CLASS][BOND],          sb)
                s3.write(row + now, col + 4, g.node[s][NODE_KEY_EDGE_CLASS][LOCAL_BRIDGE],  sb)
                s3.write(row + now, col + 5, g.node[s][NODE_KEY_EDGE_CLASS][GLOBAL_BRIDGE], sb)
                now += 1

            book.save(root + '_' + str(compNo) + '_result.xls')

            debugmsg('draw the analysis results of the target network...')
            # phase 1: draw the target network
            pos = nx.spring_layout(g, pos = nx.circular_layout(g))
            if show_the_major_result == True:
                fig_no   = 10
                bb_width = [g[s][t][EDGE_KEY_WIDTH] for (s, t) in g.edges_iter(data = False)]
                bb_color = [g[s][t][EDGE_KEY_COLOR] for (s, t) in g.edges_iter(data = False)]
                net_fig1 = pylab.figure(fig_no, figsize = (12, 8), facecolor = 'white')
                pylab.title('target network = ' + tail)
                pylab.xticks(())
                pylab.yticks(())
                net_fig1.set_tight_layout(True)
                nx.draw_networkx(g, pos = pos, linewidths = 0, width = bb_width, node_size = ns, node_color = nc, font_size = 8, edge_color = bb_color)
                pylab.savefig(root + '_' + str(compNo) + '_result.png')

            # phase 1.1: draw the detail analysis result of the target network
            if show_the_detailed_result == True:
                for i in xrange(layers):
                    l = i + 1
                    sub_edge_label = dict()
                    for s, t in g.edges_iter(data = False): sub_edge_label[(s, t)] = round(g[s][t][-l], 3)		
                    net_sub_fig = pylab.figure(fig_no + l, figsize = (12, 8), facecolor = 'white')
                    pylab.title('target network = ' + tail + ' (layer ' + str(l) + ', R1 = ' + str(round(g.graph[GRAPH_KEY_THRESHOLD_R1 + str(i + 1)], 4)) +
                                                                                   ', R2 = ' + str(round(g.graph[GRAPH_KEY_THRESHOLD_R2 + str(i + 1)], 4)) + ')')
                    pylab.xticks(())
                    pylab.yticks(())
                    net_sub_fig.set_tight_layout(True)
                    nx.draw_networkx(g, pos = pos, linewidths = 0, width = bb_width, node_size = ns, node_color = nc, font_size = 8, edge_color = bb_color)
                    nx.draw_networkx_edge_labels(g, pos = pos, edge_labels = sub_edge_label, font_size = 6)
                    pylab.savefig(root + '_' + str(compNo) + '_result_layer_' + str(l) + '.png')

            # phase 2: show betweenness centrality for edges
            if show_the_betweenness_result == True:
                eb = nx.edge_betweenness_centrality(g)
                for s, t in eb: eb[(s, t)] = round(eb[(s, t)], 3)
                bn_width = [0.5 + ((eb[(s, t)] - min(eb.values())) / scipy.std(eb.values())) for (s, t) in g.edges_iter(data = False)]            
                net_fig2 = pylab.figure(20, figsize = (12, 8), facecolor = 'white')
                pylab.title('Target network = ' + tail + ' (betweenness centrality for edges)')
                pylab.xticks(())
                pylab.yticks(())
                net_fig2.set_tight_layout(True)
                nx.draw_networkx(g, pos = pos, linewidths = 0, width = bn_width, node_size = ns, node_color = nc, font_size = 8)
                nx.draw_networkx_edge_labels(g, pos = pos, edge_labels = eb, font_size = 6)
                pylab.savefig(root + '_' + str(compNo) + '_result (edge betweenness).png')
    
            # phone 3: show pagerank-based weighting for edges
            if show_the_pagerank_result == True:
                pg = nx.Graph()
                pg.add_nodes_from(g.edges())
                for pair in pg.nodes():
                    for vertex in pair:
                        for node in g.neighbors(vertex):
                            if (vertex, node) in g.edges():
                                if not pair == (vertex, node):
                                    pg.add_edge(pair, (vertex, node))
                            else:
                                if not pair == (node, vertex):
                                    pg.add_edge(pair, (node, vertex))
                pr = nx.pagerank(pg, max_iter = 2000)
                for s, t in pr: pr[(s, t)] = round(pr[(s, t)], 4)
                pg_width = [(pr[(s, t)] - min(pr.values())) / scipy.std(pr.values()) for (s, t) in g.edges_iter(data = False)]
                net_fig3 = pylab.figure(30, figsize = (12, 8), facecolor = 'white')
                pylab.title('Target network = ' + tail + ' (pagerank-based weighting for edges)')
                pylab.xticks(())
                pylab.yticks(())
                net_fig3.set_tight_layout(True)
                nx.draw_networkx(g, pos = pos, linewidths = 0, width = pg_width, node_size = ns, node_color = nc, font_size = 8)
                nx.draw_networkx_edge_labels(g, pos = pos, edge_labels = pr, font_size = 6)
                pylab.savefig(root + '_' + str(compNo) + '_result (pagerank-based).png')

            # phase 4: show the result of network clustering
            if show_the_major_clustering_result == True:
                fig_no      = 40
                sg          = network_clustering(g, layers)
                ncc_map     = {}
                color_count = 1
                for v in g.nodes_iter(data = False):
                    if not g.node[v][NODE_KEY_GROUP_NUMBER] in ncc_map:
                        ncc_map[g.node[v][NODE_KEY_GROUP_NUMBER]] = color_count
                        color_count += 1
                ncc = [ncc_map[g.node[v][NODE_KEY_GROUP_NUMBER]] for v in g.nodes_iter(data = False)]
                net_fig4 = pylab.figure(fig_no, figsize = (12, 8), facecolor = 'white')
                pylab.title('Target network = ' + tail + ' (clustering result)')
                pylab.xticks(())
                pylab.yticks(())
                net_fig4.set_tight_layout(True)
                nx.draw_networkx(g, pos = pos, linewidths = 0, width = bb_width, node_color = ncc, vmin = min(ncc), vmax = max(ncc), cmap = pylab.cm.Dark2, font_size = 8, edge_color = bb_color)
                pylab.savefig(root + '_' + str(compNo) + '_result (network clustering).png')
                if show_the_detailed_clustering_result == True:
                    for key, sub_g in sg.items():
                        if type(sub_g) == list:
                            show_g = nx.Graph()
                            for sub_c in sub_g:
                                show_g.add_nodes_from(sub_c.nodes())
                                show_g.add_edges_from(sub_c.edges())
                        else:
                            show_g = sub_g
                        fig_no += 1
                        net_sub_fig = pylab.figure(fig_no, figsize = (12, 8), facecolor = 'white')
                        pylab.title('Target network = ' + tail + ' (' + key + ')')
                        pylab.xticks(())
                        pylab.yticks(())
                        net_sub_fig.set_tight_layout(True)
                        nx.draw_networkx(show_g, pos = pos, linewidths = 0, font_size = 8)
                        pylab.savefig(root + '_' + str(compNo) + '_result (' + key + ').png')

            pylab.show()

        return  0
        # return (g, rgs)
    else:
        return -1


def suite_experiment (suite = 'WS_SWN', onlyShow = True):    
    full_dataset = {'NWS_SWN':['nws_swn_0.0.net',    'nws_swn_0.001.net', 'nws_swn_0.002.net', 'nws_swn_0.004.net',
                               'nws_swn_0.008.net',  'nws_swn_0.016.net', 'nws_swn_0.032.net', 'nws_swn_0.064.net',
                               'nws_swn_0.128.net',  'nws_swn_0.256.net', 'nws_swn_0.384.net', 'nws_swn_0.512.net',
                               'nws_swn_0.640.net',  'nws_swn_0.768.net', 'nws_swn_0.896.net', 'nws_swn_1.0.net'],
                    
                     'WS_SWN':['ws_swn_0.0.net',     'ws_swn_0.001.net',  'ws_swn_0.002.net',  'ws_swn_0.004.net',
                               'ws_swn_0.008.net',   'ws_swn_0.016.net',  'ws_swn_0.032.net',  'ws_swn_0.064.net',
                               'ws_swn_0.128.net',   'ws_swn_0.256.net',  'ws_swn_0.384.net',  'ws_swn_0.512.net',
                               'ws_swn_0.640.net',   'ws_swn_0.768.net',  'ws_swn_0.896.net',  'ws_swn_1.0.net'],

                  'CLASSICAL':['2d.net',             'ba_sfn.net'],                
                    
                     'SOCIAL':['14p.net',            'celegans.net',      'dolphins.net',
                               'families.net',       'football.net',      'jazz.net',
                               'k-core.net',         'karate.net',        'leader.net',
                               'lesmis.net',         'prisonInter.net',   'Ragusa16.net',
                               's208.net',           'USA.net',           'women.net']}

    if onlyShow == False:
        for data_file in full_dataset[suite]:
            print "Processing " + data_file + "..."
            pylab.ion()
            experiment(data_file)
            pylab.ioff()

    nfp_shelve_path = 'network_fingerprints.obj'
    if os.path.exists(nfp_shelve_path) & os.path.isfile(nfp_shelve_path):
        sf = shelve.open(nfp_shelve_path)
        finger_prints = sf['finger_prints']
        corr_table    = sf['corr_table']
        sf.close()

        bar    = {SINK:[], GLOBAL_BRIDGE:[], LOCAL_BRIDGE:[], BOND:[]}
        labels = []
        for net_name in full_dataset[suite]:
            net_series = finger_prints[net_name[:-4] + '_1']
            labels.append(net_name[:-4])
            bar[BOND].append(net_series[0])
            bar[LOCAL_BRIDGE].append(net_series[1])
            bar[GLOBAL_BRIDGE].append(net_series[2])
            bar[SINK].append(net_series[3])
                    
        index              = numpy.arange(len(labels))
        bar[BOND]          = numpy.array(bar[BOND])
        bar[LOCAL_BRIDGE]  = numpy.array(bar[LOCAL_BRIDGE])
        bar[GLOBAL_BRIDGE] = numpy.array(bar[GLOBAL_BRIDGE])
        bar[SINK]          = numpy.array(bar[SINK])
        width              = 0.5

        net_fig1 = pylab.figure(10, figsize = (12, 8), facecolor = 'white')
        p1       = pylab.bar(index, bar[BOND],          width, color = BOND_COLOR,          edgecolor = BOND_COLOR)
        p2       = pylab.bar(index, bar[LOCAL_BRIDGE],  width, color = LOCAL_BRIDGE_COLOR,  edgecolor = LOCAL_BRIDGE_COLOR,  bottom = bar[BOND])
        p3       = pylab.bar(index, bar[GLOBAL_BRIDGE], width, color = GLOBAL_BRIDGE_COLOR, edgecolor = GLOBAL_BRIDGE_COLOR, bottom = bar[BOND] + bar[LOCAL_BRIDGE])
        p4       = pylab.bar(index, bar[SINK],          width, color = SINK_COLOR,          edgecolor = SINK_COLOR,          bottom = bar[BOND] + bar[LOCAL_BRIDGE] + bar[GLOBAL_BRIDGE])
        pylab.gca().xaxis.tick_top()                
        pylab.xticks(index + width / 2., labels, rotation = 90)
        pylab.ylabel('Percentage')
        pylab.yticks(numpy.arange(0, 1.01, 0.1))
        pylab.ylim(0, 1.)
        pylab.legend((p1[0], p2[0], p3[0], p4[0]), (BOND, LOCAL_BRIDGE, GLOBAL_BRIDGE, SINK), loc = 'lower center', fancybox = True, shadow = True, ncol = 4)
        net_fig1.set_tight_layout(True)
        for t in pylab.gca().xaxis.get_major_ticks(): t.tick1On = t.tick2On = False
        for t in pylab.gca().yaxis.get_major_ticks(): t.tick1On = t.tick2On = False                
        pylab.savefig('fingerprints_' + suite + '.png')
 
        sub_fig2, ax = pylab.subplots()        
        corr_matrix  = [[corr_table[net_name1 + '_1'][net_name2 + '_1'] for net_name2 in labels] for net_name1 in labels]
        corr_cluster = hc.dendrogram(hc.linkage(corr_matrix, method = 'centroid'), no_plot = True)
        corr_index   = corr_cluster['leaves']
        corr_result  = scipy.zeros([len(labels), len(labels)])
        for i in xrange(0, len(labels)):
            for j in xrange(0, len(labels)):
                corr_result[i, j] = corr_matrix[i][j]
        corr_result  = corr_result[corr_index, :]
        corr_result  = corr_result[:, corr_index]
        corr_labels  = [labels[i] for i in corr_index]
        corr_coef    = numpy.array(corr_result)                
        ccmap        = ax.pcolor(corr_coef, vmin = -1.0, vmax = 1.0, cmap = pylab.cm.RdBu, alpha = 0.8)
        cbar         = pylab.colorbar(ccmap)
        sub_fig2     = pylab.gcf()
        sub_fig2.set_size_inches(11, 11)
        ax.set_frame_on(False)
        ax.set_xticks(numpy.arange(corr_coef.shape[0]) + 0.5, minor = False)
        ax.set_yticks(numpy.arange(corr_coef.shape[1]) + 0.5, minor = False)
        ax.invert_yaxis()
        ax.xaxis.tick_top()
        ax.set_xticklabels(corr_labels, minor = False)
        ax.set_yticklabels(corr_labels, minor = False)
        ax.grid(False)
        pylab.xticks(rotation = 90)
        sub_fig2.set_tight_layout(True)
        for t in pylab.gca().xaxis.get_major_ticks(): t.tick1On = t.tick2On = False
        for t in pylab.gca().yaxis.get_major_ticks(): t.tick1On = t.tick2On = False
        pylab.savefig('correlation coeficient of complex networks_' + suite + '.png')

        net_fig3 = pylab.figure(12, figsize = (12, 8), facecolor = 'white')
        corr_cluster = hc.dendrogram(hc.linkage(corr_matrix, method = 'centroid'), no_plot = False, labels = labels)
        pylab.xticks(rotation = 90)
        pylab.yticks(())
        net_fig3.set_tight_layout(True)
        pylab.savefig('hierarchy of complex networks_' + suite + '.png')

        pylab.show()


def experiment(network_path, number_of_random_networks = 100, debug_msg = False):
    global path, times, debug
    
    if os.path.exists(network_path) & os.path.isfile(network_path):
        path  = network_path
        times = number_of_random_networks
        debug = debug_msg
        link_analysis()


def demo(): experiment('14p.net', debug_msg = True)
    

def main(argv = sys.argv):
    """
    主函數：剖析命令列參數。
    
    -d 表示打開除錯模式，在執行過程中顯示目前主要的執行步驟；
    -i 表示檔案及其路徑；
    -t 表示對應比較的隨機網絡個數；
    -q 表示在計算連結的強弱特性的過程中，連結的兩個端點往外擴充的層數。    
    """
    global path, times, debug, quick, separation

    warnings.filterwarnings("ignore")
    # 分析命令列參數
    try:
        opts, args = getopt.getopt(argv[1:], 'q:i:dt:')        
    except getopt.GetoptError, err:        
        print str(err)

    for o, a in opts:        
        if   o == '-d': debug = True
        elif o == '-i': path  = a
        elif o == '-t': times = int(a)
        elif o == '-q': quick, separation = (True, int(a))
        else: assert False, 'unhandled option'

    if path != None: link_analysis()


if __name__ == "__main__":
    main()
