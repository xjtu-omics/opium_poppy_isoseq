def getMaxRefEn(pairList,id):
    ret = 0
    for mp in pairList:
        ret = max(ret,mp[id])
    return ret

def getMaxIncreaseSeq(pairList):
    ena = getMaxRefEn(pairList,0)
    enb = getMaxRefEn(pairList,1)
    dp = []
    ok = []
    dp = [ [0]*(enb+1) for i in range(ena+1) ]
    ok = [ [0]*(enb+1) for i in range(ena+1) ]
    jilu = [ [0]*(enb+1) for i in range(ena+1) ]
    for mp in pairList:
        ok[mp[0]][mp[1]] = 1
    for i in range(1,ena+1):
        for j in range(1,enb+1):
            if dp[i-1][j-1] + ok[i][j] < dp[i-1][j]:
                jilu[i][j] = 1
            else:
                jilu[i][j] = 2
            dp[i][j] = max(dp[i-1][j-1] + ok[i][j], dp[i-1][j])
            if dp[i][j] < dp[i][j-1]:
                jilu[i][j] = 3
            dp[i][j] = max(dp[i][j],dp[i][j-1])
    sti,stj = ena,enb
    rel = []
    while sti>0 and stj >0:
        psti = sti
        pstj =stj
        if jilu[sti][stj] == 1:
            sti = sti - 1
        elif jilu[sti][stj] == 2:
            stj = stj - 1
            sti = sti - 1
        else:
            stj = stj - 1
        if dp[sti][stj] != dp[psti][pstj]:
            rel.append((psti,pstj))

    return dp[ena][enb],rel

def getInvSeq(pairList):
    ena = getMaxRefEn(pairList,0)
    enb = getMaxRefEn(pairList,1)
    for i in range(len(pairList)):
        pairList[i] = (pairList[i][0] , enb + 1 - pairList[i][1])
    dp = [ [0]*(enb+1) for i in range(ena+1) ]
    ok = [ [0]*(enb+1) for i in range(ena+1) ]
    jilu = [ [0]*(enb+1) for i in range(ena+1) ]
    for mp in pairList:
        ok[mp[0]][mp[1]] = 1
    for i in range(1,ena+1):
        for j in range(1,enb+1):
            if dp[i-1][j-1] + ok[i][j] < dp[i-1][j]:
                jilu[i][j] = 1
            else:
                jilu[i][j] = 2
            dp[i][j] = max(dp[i-1][j-1] + ok[i][j], dp[i-1][j])
            if dp[i][j] < dp[i][j-1]:
                jilu[i][j] = 3
            dp[i][j] = max(dp[i][j],dp[i][j-1])
    sti,stj = ena,enb
    rel = []
    while sti>0 and stj >0:
        psti = sti
        pstj = stj
        if jilu[sti][stj] == 1:
            sti = sti - 1
        elif jilu[sti][stj] == 2:
            stj = stj - 1
            sti = sti - 1
        else:
            stj = stj - 1
        if dp[sti][stj] != dp[psti][pstj]:
            rel.append((psti,pstj))
    rel.sort()
    filtRel = []
    # ====================================
    # if len(rel) > 1:
    #     for i in range(len(rel) - 1):
    #         mp1 = rel[i]
    #         mp2 = rel[i + 1]
    #         if (not ok[ mp1[0] ][ mp2[1] ]) and (not ok[ mp2[0] ][ mp1[1] ]):
    #             filtRel.append(mp1)
    #             if i == len(rel) - 2:
    #                 filtRel.append(mp2)
    if len(rel) == 2:
        mp1 = rel[0]
        mp2 = rel[1]
        if (not ok[mp1[0]][mp2[1]]) and (not ok[mp2[0]][mp1[1]]):
            filtRel = rel
    rel = filtRel
    for i in range(len(rel)):
        rel[i] = (rel[i][0] , enb + 1 - rel[i][1])
    return len(filtRel),filtRel

def getDupSeq(pairList):
    ena = getMaxRefEn(pairList,0)
    dp = [0]*(ena+1)
    rel = []
    for i in range(len(pairList)):
        mp = pairList[i]
        if dp[mp[0]] == 1:
            continue
        tmpList = []
        tmpList.append(mp)
        dp[mp[0]] = 1
        j = i+1
        while j < len(pairList):
            mpj = pairList[j]
            if mpj[0] == mp[0]:
                tmpList.append(mpj)
            j = j+1
        if len(tmpList) > 1:
            for mpj in tmpList:
                rel.append(mpj)
    return len(rel),rel