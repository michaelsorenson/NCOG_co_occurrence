library(circlize)
v4_data = read.csv('/Users/michaelsorenson/Desktop/allenlab/circlize/Verrucomicrobiota/NCOG_network_v4_circlize_edgelist_verrucomicrobiota.csv', row.names=1)
v4_positive_matr = v4_data[-c(4)]
v4_negative_matr = v4_data[-c(3)]
v4_positive_matr$Node1 = substr(v4_positive_matr$Node1, 1, 5)
v4_positive_matr$Node2 = substr(v4_positive_matr$Node2, 1, 5)
v4_negative_matr$Node1 = substr(v4_negative_matr$Node1, 1, 5)
v4_negative_matr$Node2 = substr(v4_negative_matr$Node2, 1, 5)
# chordDiagram(v4_positive_matr)
# chordDiagram(v4_negative_matr)

v9_data = read.csv('/Users/michaelsorenson/Desktop/allenlab/circlize/Verrucomicrobiota/NCOG_network_v9_circlize_edgelist_verrucomicrobiota.csv', row.names=1)
v9_positive_matr = v9_data[-c(4)]
v9_negative_matr = v9_data[-c(3)]
v9_positive_matr$Node1 = substr(v9_positive_matr$Node1, 1, 5)
v9_positive_matr$Node2 = substr(v9_positive_matr$Node2, 1, 5)
v9_negative_matr$Node1 = substr(v9_negative_matr$Node1, 1, 5)
v9_negative_matr$Node2 = substr(v9_negative_matr$Node2, 1, 5)
# chordDiagram(v9_positive_matr)
# chordDiagram(v9_negative_matr)