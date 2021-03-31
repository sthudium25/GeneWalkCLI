import pandas as pd

cbp_down = pd.read_table("/Users/Sammy/PycharmProjects/GeneWalk/venv/CBPsh5_Down_korbAllis_prelim_DEgenes.txt",
                         delimiter="\t")

print(cbp_down.head())

for col in cbp_down.columns:
    print(col)

cbp_mgi = cbp_down["MGI Gene/Marker ID"]

print(cbp_mgi.head())


cbp_filtered = cbp_down[cbp_down["MGI Gene/Marker ID"].str.contains('MGI')]
cbp_filtered = cbp_filtered["MGI Gene/Marker ID"]
print(cbp_filtered.head())
print(cbp_filtered.shape)

with open("/Users/Sammy/PycharmProjects/GeneWalk/venv/CBPsh5_Down_MGI.txt",
          mode="w") as f:
    f.write(cbp_filtered)
