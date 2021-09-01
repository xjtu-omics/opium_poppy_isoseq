import pandas as pd
from collections import defaultdict
#返回一个适用于seaborn的barplot的累积百分比条形图
def countDf2PercentageDf(data,xName,colorName,countName):
#因为seaborn的barplot没有直接将数值累加的API,所以只能将原始数据进行转换
#先将原始数据转换成百分比数据以后,按照字典序的顺序反向累加,因为seaborn采取
#从前往后覆盖上色的的策略,所以要把最后上色的放到累加的第一顺序
    xSet = data.loc[:,[xName]].squeeze().drop_duplicates().tolist()
    colorSet = data.loc[:,[colorName]].squeeze().drop_duplicates().tolist()
    xCount = defaultdict(int)
    for i in range(len(xSet)):
        x = xSet[i]
        rowIndex = (data[[xName]]==x).squeeze()
        xCount[x] = data.loc[ rowIndex,[countName] ].sum().at[countName]
    percCol =  []
    for xi in range(data.shape[0]):
        tmpCount = 0
        xiX = data.at[xi,xName]
        xiColor = data.at[xi,colorName]
        xiXIndex = (data[[xName]]==xiX).squeeze()
        xiXSubData = data.loc[xiXIndex,:]
        xiXSubData.index = xiXSubData[[colorName]].squeeze()
        for i in range(len(colorSet)):
            j = len(colorSet)-i-1
            tmpCount = tmpCount + xiXSubData.at[colorSet[j],countName]
            if colorSet[j] == xiColor:
                break
        tmpXCount = xCount[ data.at[ xi,xName ] ]
        percCol.append(tmpCount/tmpXCount)
    percCol = pd.Series(percCol)
    relpd = pd.concat([data,pd.Series(percCol)], axis=1)
    colNames = data.columns.values.tolist()
    colNames.append('perc')
    relpd.columns = colNames
    return(relpd)
