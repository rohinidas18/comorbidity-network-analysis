import pandas as pd

xls = pd.ExcelFile('net5 - ra.xlsx')
df1 = pd.read_excel(xls, 'dis', header=None)
#df2 = pd.read_excel(xls, 'common', header=None)
df3 = pd.read_excel(xls, 'tb', header=None)

tw = pd.read_csv('string_interactions_short - ra_default node.csv')

def conditions(s):
    if s in list(df1[0]):
        return "dis"
    #elif s in list(df2[0]):
        #return "common"
    elif s in list(df3[0]):
        return "tb"

tw['type'] = tw['name'].apply(lambda row: conditions(row))
tw.to_csv('string_interactions_short - ra_c.csv', index=False)
