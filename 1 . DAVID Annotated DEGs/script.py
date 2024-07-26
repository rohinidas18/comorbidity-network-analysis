import pandas as pd

# f = open("siglist56.txt", "r")
# fcont = f.readlines()
# genenames = {}

# for i in fcont:
#     ty = i.split("\t")
#     genenames[ty[0]] = ty[1]

# df = pd.DataFrame(data=genenames, index=[0])
# df = df.transpose()
# df.to_excel('Gene ID vs Name - 5.xlsx')




f = open("siglist56.txt", "r")
fcont = f.readlines()
genenames = []

for i in fcont:
    ty = i.split("\t")
    genenames.append(ty)



#df = pd.DataFrame(data=genenames, index=[0])
#df = df.transpose()
# df.to_excel('Gene ID vs Name - 5.xlsx')