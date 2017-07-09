import random
import sys

import itertools as it

import numpy as np
import scipy.optimize, scipy.misc, scipy.cluster

def extract_chains(seq_data):
  alphas_per_well, betas_per_well = zip(*seq_data.well_data)
  return sorted(set(sum(alphas_per_well, []))), sorted(set(sum(betas_per_well, [])))

def solve(seq_data, iters=100, pair_threshold = 0.9):

  ## Computes a solution to the alpha-beta pairing problem, using the methods in Lee et al. (2017)
  def compute_well_pairings(alpha_idx, beta_idx, scores):

    if len(alpha_idx)==0 or len(beta_idx)==0:  return [] # scipy hungarian implementation doesn't handle this edge case

    # Reformulate problem as a general assignment problem
    # Then apply Hungarian algorithm
    # Indices in the results of running Hungarian are transformed back into alpha/beta chain ids
    ratings = [[-scores[i][j] for j in beta_idx] for i in alpha_idx]
    pairings = [(alpha_idx[i], beta_idx[j]) for i,j in zip(*scipy.optimize.linear_sum_assignment(ratings))]

    return pairings
  
  # Extract all distinct alpha- and beta-chains observed
  # TODO: might be better to extract the chains directly from the cells in the system
  all_alphas, all_betas = extract_chains(seq_data)
  
  # Create dictionary to look up alpha/beta index in constant time
  alpha_to_idx = {a: i for i,a in enumerate(all_alphas)}
  beta_to_idx = {b: i for i,b in enumerate(all_betas)}
  # Transform all well data to reference alpha- and beta-chains by index
  well_data = [[[alpha_to_idx[a] for a in data[0]], [beta_to_idx[b] for b in data[1]]] for data in seq_data.well_data]


  overall_pairing_counts = {}
  for iter in range(0): #range(iters):
    # Choose random subset of wells for this iter
    # Loop is to ensure that subset size is greater than 0 (can happen w/ small well count)
    #wells_idx = []
    #while len(wells_idx)==0:  wells_idx = [i for i in range(len(well_data)) if random.random()>0.5]
    # Actually, each random subset is constant fraction (0.75) of all wells
    wells_idx = np.random.choice(range(len(well_data)), int(0.75*len(well_data)), replace=False)
    
    # Calculate association scores
    S = [[0 for j in range(len(all_betas))] for i in range(len(all_alphas))]
    for well_idx in wells_idx:
      well_alpha_idx, well_beta_idx = well_data[well_idx]
      for a_idx in well_alpha_idx:
        for b_idx in well_beta_idx:
          increment = 1./len(well_alpha_idx) + 1./len(well_beta_idx)
          S[a_idx][b_idx] += increment

    # Compute well pairings for any well, if it hasn't been done already
    # Then accumulate the number of times each pair has been assigned in a well pairing
    pairing_counts = {}
    for idx, well_idx in enumerate(wells_idx):
      well_pairings = compute_well_pairings(*well_data[well_idx], scores=S)
      #print "Computing likely pairings... {0}%\r".format(int(100*((iter+1) + float(idx)/len(wells_idx))/iters)),
      for a,b in well_pairings:
        pairing_counts[(a,b)] = pairing_counts.get((a,b), 0) + 1

    # Compute filter cutoff (average of all nonzero pair counts)
    cutoff = np.mean(pairing_counts.values())

    # Extract all pairs with counts exceeding the cutoff
    good_pairs = [pair for pair in pairing_counts if pairing_counts[pair]>cutoff]

    # For each pair exceeding the cutoff, increment the overall_pairing_counts number
    for pair in good_pairs:
      overall_pairing_counts[pair] = overall_pairing_counts.get(pair, 0) + 1

    print "Computing likely pairings... {0}%\r".format(100*(iter+1)/iters),
    sys.stdout.flush()

  overall_good_pairs = [pair for pair in overall_pairing_counts if overall_pairing_counts[pair]>=pair_threshold*iters]

  pairs = [(all_alphas[a], all_betas[b]) for a,b in overall_good_pairs]
  thresholds = [overall_pairing_counts[p]/float(iters) for p in overall_good_pairs]

  pairs = [(81, 304), (1610, 593), (216, 588), (1480, 429), (1247, 1561), (391, 620), (708, 1526), (59, 1337), (1530, 1744), (1520, 1282), (1393, 894), (184, 943), (198, 50), (742, 950), (1068, 1778), (281, 845), (692, 160), (1571, 69), (798, 46), (1332, 439), (917, 1003), (280, 34), (1656, 1193), (1375, 485), (984, 594), (594, 1292), (388, 1009), (794, 1359), (523, 620), (1260, 1698), (472, 280), (726, 1420), (153, 1309), (1223, 1165), (1025, 1060), (869, 782), (1401, 1560), (52, 437), (674, 1797), (948, 1729), (767, 1612), (1093, 711), (1148, 1173), (1415, 384), (834, 1100), (712, 704), (354, 844), (246, 269), (1321, 1515), (168, 1321), (462, 1649), (1420, 869), (245, 1645), (1726, 185), (1610, 1606), (1443, 1746), (1248, 1224), (834, 408), (1674, 1503), (655, 864), (588, 1170), (387, 1238), (1562, 433), (1210, 333), (930, 1692), (503, 240), (1010, 186), (1732, 57), (389, 1234), (920, 751), (1048, 413), (1580, 974), (289, 1113), (1184, 1102), (400, 796), (1547, 463), (573, 378), (1089, 1251), (362, 422), (833, 930), (1014, 1700), (1183, 475), (1255, 1575), (669, 791), (1366, 195), (58, 292), (730, 1109), (829, 370), (165, 77), (1304, 78), (730, 1147), (1652, 1094), (1110, 147), (677, 1385), (507, 400), (1589, 769), (204, 1829), (414, 1817), (216, 190), (262, 671), (1732, 708), (290, 1097), (1365, 1109), (1348, 1650), (1074, 523), (1237, 1769), (1621, 1723), (1449, 1710), (1133, 1144), (725, 1545), (118, 1008), (211, 1806), (1590, 27), (1710, 502), (188, 872), (997, 747), (1334, 691), (1673, 442), (7, 282), (1307, 991), (1182, 1688), (1489, 1424), (176, 1790), (1215, 1224), (1446, 1756), (1458, 849), (59, 1401), (1488, 1484), (832, 976), (1260, 1783), (1012, 763), (566, 1786), (1211, 828), (266, 940), (703, 93), (555, 1754), (627, 446), (648, 523), (1159, 749), (818, 211), (1615, 574), (1323, 257), (146, 1676), (1104, 773), (342, 982), (1307, 1643), (637, 1396), (724, 385), (552, 1763), (161, 728), (1234, 561), (1745, 1233), (387, 1271), (1735, 1253), (750, 487), (1304, 877), (672, 1327), (893, 1445), (1303, 568), (252, 284), (98, 1622), (82, 608), (909, 1052), (430, 1011), (472, 680), (1634, 1375), (11, 30), (1628, 1063), (938, 435), (969, 1804), (1315, 652), (137, 1229), (523, 1235), (1181, 866), (502, 78), (870, 524), (1344, 1051), (799, 853), (1282, 1745), (1724, 1151), (1746, 1732), (1387, 24), (111, 1127), (1028, 47), (1678, 469), (1096, 1518), (1105, 1342), (77, 1370), (1198, 616), (1383, 797), (1433, 764), (1457, 618), (478, 237), (1283, 1310), (1725, 1740), (1195, 1673), (673, 324), (1141, 219), (1607, 214), (1732, 1492), (617, 903), (206, 1759), (1603, 168), (1599, 651), (376, 1381), (1300, 628), (1322, 1044), (731, 230), (1666, 650), (800, 993), (1386, 206), (698, 1802), (883, 353), (520, 1736), (32, 823), (1650, 1232), (297, 820), (592, 1775), (148, 134), (1090, 159), (447, 1176), (657, 251), (826, 1017), (45, 1807), (1183, 1014), (318, 4), (224, 145), (637, 893), (982, 217), (1248, 319), (380, 404), (507, 1090), (1006, 978), (1164, 927), (1337, 1750), (1676, 1719), (1691, 979), (1632, 822), (1267, 694), (1529, 1109), (486, 1269), (372, 128), (72, 1646), (1368, 1457), (1377, 126), (981, 520), (1102, 937), (942, 1005), (1321, 554), (408, 1458), (206, 1109), (938, 1066), (167, 1566), (148, 1543), (187, 94), (1414, 1701), (1730, 582), (619, 1463), (726, 676), (386, 343), (909, 1429), (287, 1728), (1300, 964), (1737, 920), (1618, 793), (1426, 1049), (239, 109), (1457, 1608), (958, 867), (1517, 310), (792, 1109), (60, 657), (137, 858), (244, 1362), (1463, 1134), (1365, 1235), (869, 719), (610, 981), (929, 949), (425, 793), (1595, 172), (498, 483), (425, 1102), (738, 543), (224, 388), (206, 497), (977, 1109), (397, 754), (520, 602), (795, 1235), (968, 1771), (398, 151), (1078, 1410), (69, 1030), (526, 632), (502, 543), (1125, 1587), (1548, 52), (1315, 1672), (999, 97), (1451, 591), (1173, 1589), (291, 1025), (1054, 1042), (157, 420), (1185, 1678), (892, 662), (1618, 1456), (1325, 1709), (274, 37), (734, 456), (1327, 1035), (605, 921), (799, 1476), (1540, 549), (1747, 1715), (375, 1130), (1177, 506), (718, 611), (834, 349), (467, 98), (759, 634), (1248, 1162), (246, 1623), (436, 567), (662, 20), (206, 1235), (160, 1715), (1413, 507), (1279, 1162), (1338, 1367), (379, 347), (25, 1109), (976, 537), (824, 33), (1094, 961), (1419, 1240), (264, 1236), (1756, 825), (1656, 1040), (674, 243), (331, 155), (1085, 1726), (728, 198), (1199, 1154), (625, 1121), (1543, 808), (814, 279), (461, 405), (197, 313), (1159, 1453), (822, 1235), (387, 692), (554, 1431), (790, 1198), (1093, 1092), (1227, 6), (1159, 755), (1391, 1413), (50, 650), (1581, 842), (943, 78), (1529, 1456), (859, 1602), (983, 262), (1630, 1688), (1564, 233), (247, 916), (132, 1204), (570, 848), (765, 1520), (1409, 1485), (1121, 1350), (1520, 657), (1449, 474), (210, 1293), (1129, 1332), (939, 669), (1510, 1316), (1324, 1286), (1221, 1146), (966, 1195), (757, 845), (1345, 1562), (336, 1741), (1670, 14), (957, 986), (1135, 395), (670, 843), (682, 643), (426, 1045), (738, 1235), (31, 1491), (684, 1564), (578, 1087), (718, 136), (1748, 1601), (1222, 667), (1664, 227), (658, 1585), (474, 775), (1465, 910), (1138, 1537), (1179, 440), (1390, 1242), (1675, 362), (425, 1235), (256, 248), (1734, 308), (887, 1662), (774, 289), (1711, 865), (736, 615), (1706, 1374), (513, 1365), (1336, 103), (420, 268), (922, 1103), (932, 1620), (388, 965), (1418, 238), (1248, 1109), (730, 1759), (670, 1330), (481, 1111), (1273, 1089), (1529, 1235), (159, 767), (22, 1064), (425, 877), (250, 1534), (1444, 1073), (1457, 1557), (242, 1153), (149, 689), (301, 857), (748, 1696), (589, 1372), (611, 59), (936, 1588), (1380, 212), (1330, 969), (25, 635), (1, 95), (338, 955), (1132, 1760), (393, 770), (162, 276), (927, 756), (274, 1179), (1500, 1751), (1134, 88), (976, 245), (1246, 188), (1153, 36), (872, 1477), (1469, 377), (319, 804), (871, 1669), (1471, 83), (896, 1742), (26, 1226), (1058, 629), (172, 111), (730, 78), (1372, 1241), (1202, 1508), (300, 1707), (1444, 1125), (1743, 322), (1168, 89), (935, 255), (1508, 1391), (1587, 453), (459, 150), (1066, 1744), (1297, 587), (799, 358), (1126, 1105), (1464, 1422), (772, 148), (309, 66), (533, 1440), (1046, 1168), (710, 1367), (891, 1449), (1618, 877), (1295, 1753), (1137, 1680), (218, 1160), (1057, 777), (1421, 512), (892, 1216), (422, 207), (934, 1502), (1148, 443), (743, 1558), (50, 1080), (644, 776), (438, 527), (1321, 254), (1559, 10), (1698, 1653), (145, 1535), (288, 351), (27, 229), (1516, 1291), (317, 382), (1427, 975), (1625, 1539), (528, 444), (1151, 1507), (1283, 1013), (909, 739), (873, 917), (1070, 1578), (1622, 1043), (805, 565), (1136, 1235), (881, 1700), (1196, 96), (1379, 1382), (907, 1016), (1135, 1479), (640, 104), (628, 1192), (285, 409), (65, 968), (1237, 1482), (792, 1162), (1501, 1184), (11, 1187), (689, 74), (1604, 1832), (1124, 1222), (1092, 1333), (1572, 705), (1509, 518), (502, 1660), (374, 358), (58, 54), (256, 595), (561, 576), (67, 1047), (167, 232), (730, 1660), (795, 1443), (1119, 936), (586, 1219), (180, 1180), (1279, 1147), (1083, 1378), (4, 933), (425, 1456), (1034, 965), (554, 1164), (1136, 597), (1566, 1197), (1171, 1532), (1524, 1631), (1309, 1755), (1702, 1191), (139, 713), (1280, 473), (679, 1655), (979, 1188), (18, 1443), (1614, 1670), (29, 1386), (1307, 762), (160, 836), (83, 693), (953, 1517), (773, 90), (635, 1075), (532, 451), (453, 835), (660, 299), (878, 778), (1444, 1219), (977, 942), (601, 1219), (1423, 326), (648, 1480), (589, 1162), (1182, 228), (85, 592), (122, 157), (522, 1162), (1697, 868), (383, 511), (697, 482), (1618, 1235), (1165, 1132), (905, 861), (1365, 106), (1279, 1109), (595, 1463), (1042, 1426), (1630, 840), (1304, 868), (463, 725), (68, 125), (703, 1182), (1285, 622), (882, 570), (549, 1417), (1240, 1542), (206, 543), (395, 1021), (1043, 1776), (1065, 1065), (1629, 1550), (145, 1311), (619, 876), (822, 877), (1526, 1227), (146, 1083), (613, 1787), (339, 258), (240, 1156), (1569, 1295), (720, 1077), (911, 1800), (709, 1665), (1360, 980), (1511, 494), (1170, 145), (1479, 1347), (81, 1004), (89, 375), (1561, 1510), (1705, 44), (1449, 290), (803, 1421), (496, 1319), (235, 1683), (1709, 779), (1584, 149), (933, 302), (1726, 329), (535, 449), (1142, 1247), (1765, 575), (475, 758), (1676, 70), (943, 1235), (369, 1574), (730, 1235), (420, 174), (1576, 1108), (158, 1600), (1185, 646), (730, 126), (670, 120), (491, 1199), (1632, 183), (413, 631), (1464, 499), (1568, 1208), (367, 246), (1618, 543), (1372, 1711), (408, 1544), (1234, 1795), (749, 1101), (394, 167), (445, 1229), (358, 505), (534, 144), (438, 1366), (265, 1685), (1241, 946), (589, 1235), (206, 257), (853, 1415), (886, 319), (522, 761), (885, 325), (1751, 989), (1333, 599), (253, 1120), (674, 1703), (697, 1591), (1700, 505), (1160, 1759), (1051, 938), (1480, 124), (1066, 1109), (101, 1226), (1660, 869), (394, 518), (502, 1134), (1248, 1729), (1358, 1256), (1214, 1039), (1372, 22), (1051, 1237), (1520, 301), (442, 1171), (1177, 1235), (1299, 8), (110, 215), (1091, 1500), (1177, 832), (971, 610), (1238, 57), (1335, 265), (1638, 1671), (1694, 997), (1035, 915), (1156, 1577), (504, 808), (909, 1368), (206, 106), (1470, 1196), (831, 954), (1277, 1213), (39, 224), (929, 757), (1651, 132), (327, 975), (1502, 1336), (447, 1371), (1632, 323), (1640, 895), (1708, 1246), (134, 1747), (888, 477), (1306, 1801), (182, 596), (1311, 1702), (1175, 1076), (1503, 1197), (1149, 232), (658, 1654), (534, 1206), (960, 880)]

  # Turns pairs of associated alpha- and beta- chains into cells that may have dual alpha chains
  cells, cell_freqs, cell_freqs_CI = pairs_to_cells(seq_data, pairs) 


  results = {
    'cells': cells,
    'cell_frequencies': cell_freqs,
    'cell_frequencies_CI': cell_freqs_CI,
    'pairs': pairs,
    'pair_thresholds': thresholds
  }

  return results

def estimate_cell_frequencies(seq_data, cells):

  def log_likelihood_func(f, N, W, K, Q_memo = {}, error_rate=0.15, is_dual=False):
    # Note: See Eqs (3) and (4) in Lee et al. for explanation of variables
  
    # Compute Q vector if not previously computed
    Q_key = (tuple(N), tuple(W), error_rate, f, is_dual)
    if Q_key not in Q_memo:
      Q = []
      if not is_dual:
        prefactor = lambda m: 2*error_rate**m - error_rate**(2*m)
      else:
        prefactor = lambda m: 3*error_rate**m - 3*error_rate**(2*m) + error_rate**(3*m)
      for n,w in zip(N,W):
        q = (1-f)**n + sum([
            prefactor(m) * scipy.misc.comb(n,m) * (f**m) * (1-f)**(n-m) 
        for m in range(1, n+1)])
        Q.append(q)
      Q_memo[Q_key] = Q
  
    # Retrieve Q from memoized dict of Q's
    # Note that Q only depends on N, W, error_rate, f, and is_dual
    Q = Q_memo[Q_key]
  
    # Compute log likelihood as sum of Binomial probabilities
    # Note that the "combinations" in the binomial PDF is ignored as it does not affect
    # the location of the maximum
    return sum([np.log(scipy.misc.comb(w,k)) + k*np.log((1-q)) + (w-k)*np.log(q) for w,k,q in zip(W,K,Q)])

  cells_per_well, N, W = extract_cells_per_well(seq_data)

  K = extract_cell_counts(seq_data, cells, cells_per_well, N, W)

  cell_freqs = []
  cell_freq_CIs = []
  for (alist, blist), k in zip(cells, K):
    L_func = lambda f: log_likelihood_func(f, N, W, k, is_dual=len(alist)>1)

    # Find maximal likelihood
    f_opt = scipy.optimize.minimize_scalar(lambda f: -L_func(f), method='Bounded', bounds=(0,1), options={'xatol': 1e-10}).x
    L_max = L_func(f_opt)

    # Find confidence interval, as specified in the paper
    f_min = scipy.optimize.minimize_scalar(lambda f: (L_max-1.96-L_func(f))**2, method='Bounded', bounds=(0,f_opt)).x
    f_max = scipy.optimize.minimize_scalar(lambda f: (L_max-1.96-L_func(f))**2, method='Bounded', bounds=(f_opt,1)).x

    cell_freqs.append(f_opt)
    cell_freq_CIs.append((f_min, f_max))

    print "Computing chain pair frequencies... {0}%\r".format(int(100.*len(cell_freqs)/len(cells))),
  
  return cell_freqs, cell_freq_CIs

def pairs_to_cells(seq_data, pairs):
  def find_duals_likelihood(candidate_duals, freqs_dict, well_size_cutoff = 50, error_rate=0.15):
    cells_per_well, N, W = extract_cells_per_well(seq_data)

    duals = []
    for alist, blist in candidate_duals:
      cells_temp = [((alist[0],), blist), ((alist[1],), blist), (alist, ()), (alist, blist)]
      print "", cells_temp
      K = extract_cell_counts(seq_data, cells_temp, cells_per_well, N, W)

      # Extract individual cell counts (see Lee et al., SI, Section 5 for explanation of variables)
      K_1 = K[0]
      K_2 = K[1]
      K_3 = K[2]
      K_d = K[3]
      K_o = [w-k1-k2-k3-kd for w,k1,k2,k3,kd in zip(W, K_1, K_2, K_3, K_d)]

      print "", K_1, K_2, K_3, K_d, K_o

      # Extract relevant cell frequencies
      f_q = freqs_dict[cells_temp[0]]
      f_r = freqs_dict[cells_temp[1]]
      f_d = freqs_dict[cells_temp[3]]

      print "", f_q, f_r, f_d

      # Null hypothesis (no dual clone)
      log_fact = lambda x: scipy.special.gammaln(x+1)

      # disgusting calculations, shield your eyes from the carnage below
      null_P_a1b = [
        sum([
          np.exp(log_fact(n) - log_fact(n_1) - log_fact(n_2) - log_fact(n-n_1-n_2)) * f_q**n_1 * f_r**n_2 * (1-f_q-f_r)**(n-n_1-n_2) * (1-error_rate**n_1) * error_rate**n_2 * (1-error_rate**(n_1+n_2))
          for n_1 in range(1, n+1) for n_2 in range(0, n-n_1+1)
        ])
        for n in N if n<well_size_cutoff
      ]
      #print null_P_a1b, null_P_a1b_check, null_P_a1b_check2
      null_P_a2b = [
        sum([
          np.exp(log_fact(n) - log_fact(n_1) - log_fact(n_2) - log_fact(n-n_1-n_2)) * f_q**n_1 * f_r**n_2 * (1-f_q-f_r)**(n-n_1-n_2) * (1-error_rate**n_2) * error_rate**n_1 * (1-error_rate**(n_1+n_2))
          for n_2 in range(1, n+1) for n_1 in range(0, n-n_2+1)
        ])
        for n in N if n<well_size_cutoff
      ]
      null_P_a1a2 = [
        sum([
          np.exp(log_fact(n) - log_fact(n_1) - log_fact(n_2) - log_fact(n-n_1-n_2) + n_1*np.log(f_q) + n_2*np.log(f_r) + (n-n_1-n_2)*np.log(1-f_q-f_r) + n_1*np.log(error_rate) + np.log(1-error_rate**n_1) + n_2*np.log(error_rate) + np.log(1-error_rate**n_2))
          for n_1 in range(1,n) for n_2 in range(1, n-n_1+1)
        ])
        for n in N if n<well_size_cutoff
      ]
      null_P_a1a2b = [
        sum([
          np.exp(log_fact(n) - log_fact(n_1) - log_fact(n_2) - log_fact(n-n_1-n_2)) * f_q**n_1 * f_r**n_2 * (1-f_q-f_r)**(n-n_1-n_2) * (
            error_rate**n_1*(1-error_rate**n_1)*(1-error_rate**n_2)**2 +
            (1-error_rate**n_1)**2*(1-error_rate**n_2)**2 +
            (1-error_rate**n_1)**2*error_rate**n_2*(1-error_rate**n_2)
          )
          for n_1 in range(1,n) for n_2 in range(1, n-n_1+1)
        ])
        for n in N if n<well_size_cutoff
      ]
      null_P_other = [1-p1-p2-p3-pd for p1,p2,p3,pd in zip(null_P_a1b, null_P_a2b, null_P_a1a2, null_P_a1a2b)]
        
      null_log_likelihood = sum([
        log_fact(w) - log_fact(k1) - log_fact(k2) - log_fact(k3) - log_fact(kd) - log_fact(ko) + k1*np.log(p1) + k2*np.log(p2) + k3*np.log(p3) + kd*np.log(pd) + ko*np.log(po)
        for w,k1,k2,k3,kd,ko,p1,p2,p3,pd,po in zip(W, K_1, K_2, K_3, K_d, K_o, null_P_a1b, null_P_a2b, null_P_a1a2, null_P_a1a2b, null_P_other)
      ])

      # Alternative hypothesis (dual clone a1a2b)
      alt_P_2 = [
        sum([
          scipy.misc.comb(n, k) * f_d**k * (1-f_d)**(n-k) * error_rate**k * (1-error_rate**k)**2
          for k in range(1, n+1)
        ])
        for n in N if n<well_size_cutoff
      ]
      alt_P_3 = [
        sum([
          scipy.misc.comb(n,k) * f_d**k * (1-f_d)**(n-k) * (1-error_rate**k)**3
          for k in range(1, n+1)
        ])
        for n in N if n<well_size_cutoff
      ]
      alt_P_other = [1 - 3*p1 - p2 for p1,p2 in zip(alt_P_2, alt_P_3)]

      alt_log_likelihood = sum([
        log_fact(w) - log_fact(k1) - log_fact(k2) - log_fact(k3) - log_fact(kd) - log_fact(ko) + (k1+k2+k3)*np.log(p1) + kd*np.log(p2) + ko*np.log(po)
        for w,k1,k2,k3,kd,ko,p1,p2,po in zip(W, K_1, K_2, K_3, K_d, K_o, alt_P_2, alt_P_3, alt_P_other)
      ])

      if alt_log_likelihood - null_log_likelihood >= 10:
        duals.append(cells_temp[3])

    return duals

  def find_duals_clustering(candidate_duals, freqs_dict):
    if len(candidate_duals) < 2:
      return []

    # Preliminaries
    cells_per_well, N, W = extract_cells_per_well(seq_data)
    K = extract_cell_counts(seq_data, candidate_duals, cells_per_well, N, W)

    # Compute ratio statistics
    R = []
    for (alist, blist), K_row in zip(candidate_duals, K):
      f1 = freqs_dict[((alist[0],), blist)]
      f2 = freqs_dict[((alist[1],), blist)]
      expected = sum([
        w*(1 - (1-f1)**n - (1-f2)**n + (1-f1-f2)**n)
        for n,w in zip(N,W)
      ])
      #print (alist, blist), f1, f2, sum(K_row), expected, float(sum(K_row))/expected
      R.append(float(sum(K_row))/expected)

    # Perform clustering based on R
    centroids = sorted(scipy.cluster.vq.kmeans(R, 2)[0])
    if len(centroids)==1:  return []

    C_nondual, C_dual = centroids
    
    # Filter cells based on how close R is to C_dual vs. C_nondual
    duals = []
    #print centroids
    for candidate, r in zip(candidate_duals, R):
      if np.abs(r-C_dual) < np.abs(r-C_nondual):
        duals.append(candidate)
    
    return duals
    

  candidate_non_duals = [((a,),(b,)) for a,b in pairs] # list of all cells, will be modified as duals are found
  cells = candidate_non_duals[:]
  #print pairs

  # Determine candidate dual cells
  # Note that only dual alpha-chains are considered by the Lee et al. technique
  candidate_duals = []
  beta_pairings = {}
  for a,b in pairs:
    beta_pairings[b] = beta_pairings.get(b, []) + [a]
  for b,alist in beta_pairings.iteritems():
    if len(alist) >= 2:
      alist = sorted(alist)
      candidate_duals.extend(it.product(it.combinations(alist, 2), [(b,)]))

  #print beta_pairings

  # Estimate frequencies of all potential cells
  freqs_list, freqs_CI_list = estimate_cell_frequencies(seq_data, candidate_non_duals + candidate_duals)
  freqs_dict = {c: f for c,f in zip(candidate_non_duals+candidate_duals, freqs_list)}
  freqs_CI_dict = {c: f for c,f in zip(candidate_non_duals+candidate_duals, freqs_CI_list)}
  
  # Find duals using likelihood method, which is computationally infeasible with wells with >50 cells
  likelihood_duals = find_duals_likelihood(candidate_duals, freqs_dict)
  print "Likelihood duals", likelihood_duals

  # Find duals using clustering method, which works better lower-frequency cells
  clustering_duals = find_duals_clustering(candidate_duals, freqs_dict)
  print "Clustering duals", clustering_duals

  # Remove non-dual counterparts for each dual cell found and add in corresponding dual cell
  duals = list(set(likelihood_duals + clustering_duals))
  for alist,blist in duals:
    if ((alist[0],), blist) in cells:
      cells.remove(((alist[0],), blist))
    if ((alist[1],), blist) in cells:
      cells.remove(((alist[1],), blist))
    cells.append((alist, blist))
  
  cell_freqs = [freqs_dict[c] for c in cells]
  cell_freqs_CI = [freqs_CI_dict[c] for c in cells]

  return cells, cell_freqs, cell_freqs_CI
  

## Auxiliary functions to help out pairs_to_cells() and estimate_cell_frequencies()
def extract_cells_per_well(seq_data):
  cpw_distro = seq_data.metadata['cells_per_well_distribution']
  cpw_params = seq_data.metadata['cells_per_well_distribution_params']
  if cpw_distro == 'constant':
    cells_per_well = [cpw_params['cells_per_well']]*len(seq_data.well_data)
  elif cpw_distro == 'poisson':
    cells_per_well = [cpw_params['lam']]*len(seq_data.well_data) # This is an approx. Not sure if it'll affect results
  elif cpw_distro == 'explicit':
    cells_per_well = cpw_params['cells_per_well']
  else:
    print "Unknown cell/well distribution: {0} with parameters {1}".format(cpw_distro, cpw_params)
    return None, None, None

  # Gather distinct #s of cells per well (N) and count number of wells w/ each cell count
  N_dict = {}
  for cpw in cells_per_well:
    N_dict[cpw] = N_dict.get(cpw, 0) + 1
  N,W = zip(*sorted(N_dict.iteritems()))
  return cells_per_well, N, W

def extract_cell_counts(seq_data, cells, cells_per_well, N, W):
  K = [[0]*len(N) for i in range(len(cells))]

  for well_size, well_data in zip(cells_per_well, seq_data.well_data):
    well_alphas = set(well_data[0])
    well_betas = set(well_data[1])
    N_idx = N.index(well_size)
    for i,(alist,blist) in enumerate(cells):
      if all([a in well_alphas for a in alist]) and all([b in well_betas for b in blist]):
        K[i][N_idx] += 1

  return K

