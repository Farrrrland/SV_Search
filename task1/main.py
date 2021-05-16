import time
with open("sv.fasta",'r') as sv_fa:
    SV = {}
    for line in sv_fa:
        if line.startswith(">"):
            name = line.rstrip("\n")
            name = name[1: ]
            SV[name] = ""
        else:
            SV[name] = SV[name] + line.rstrip("\n")
sv_fa.close()

with open("ref.fasta",'r') as ref_fa:
    REF = {}
    for line in ref_fa:
        if line.startswith(">"):
            name = line.rstrip("\n")
            name = name[1: ]
            REF[name] = ""
        else:
            REF[name] = REF[name] + line.rstrip("\n")
ref_fa.close()

fout = open("sv.bed", 'w')

TRA = []
sv = []
ref = []
min = 50
max = 1001

def match (i, j, len1, len2):
    cnt = 0
    while i < len1 and j < len2 and cnt < max:
        if sv[i] == ref[j] or sv[i] == 'N' or ref[j] == 'N':
            cnt += 1
            i += 1
            j += 1
        else:
            break
    return cnt


def greedy_search(i, j, len1, len2):
    case = 0
    max_match = 0
    refpos1 = -1
    refpos2 = -1
    svpos1 = -1
    svpos2 = -1
    # traverse sv
    for idx in range(min, max):
        if i+idx >= len1:
            break
        if sv[i+idx] == ref[j] or sv[i+idx] == 'N' or ref[j] == 'N':
            # match backwards
            cnt = match(i + idx, j, len1, len2)
            if cnt > max_match:
                # determine dup or ins
                max_match = cnt
                if hash(sv[i: i + idx]) == hash(sv[i - idx: i]):
                    case = 2    # duplicate
                    svpos1 = i
                    svpos2 = i + idx
                    refpos1 = j
                    refpos2 = j + idx
                    if cnt == max:
                        return case, refpos1, refpos2, svpos1, svpos2, max_match
                else:
                    case = 1    # insert
                    svpos1 = i
                    svpos2 = i + idx
                    refpos1 = j
                    refpos2 = j + idx 
                    if cnt == max:
                        return case, refpos1, refpos2, svpos1, svpos2, max_match
    # traverse ref
    for idx in range(min, max):
        if j+idx >= len2:
            break
        if ref[j+idx] == sv[i] or sv[i] == 'N' or ref[j+idx] == 'N':
            cnt = match(i, j + idx, len1, len2)
            if cnt > max_match:
                max_match = cnt
                case = 3    # delete
                refpos1 = j
                refpos2 = j + idx
                svpos1 = i
                svpos2 = i + idx
                if cnt == max:
                    return case, refpos1, refpos2, svpos1, svpos2, max_match
    # ref & sv traverse at the same time
    for idx in range(min, max):
        if j+idx >= len2 or i+idx >= len1:
            break
        if ref[j+idx] == sv[i+idx] or sv[i+idx] == 'N' or ref[j+idx] == 'N':
            cnt = match(i+idx, j+idx, len1, len2)
            if cnt > max_match:
                max_match = cnt
                is_inv = True
                k = 1
                while i + idx - k >= i:
                    if sv[i + idx - k] == 'N' or ref[j + k - 1] == 'N':
                        is_inv = True
                    elif sv[i + idx - k] == 'A' and ref[j + k - 1] != 'T' or sv[i + idx - k] == 'T' and ref[j + k - 1] != 'A' or sv[i + idx - k] == 'C' and ref[j + k - 1] != 'G' or sv[i + idx - k] == 'G' and ref[j + k - 1] != 'C':
                        is_inv = False
                        break
                    k += 1
                if is_inv == True:
                    case = 4    # inverse
                else:
                    case = 5    # TRA
                refpos1 = j
                refpos2 = j + idx
                svpos1 = i
                svpos2 = i + idx
                if cnt == max:
                    return case, refpos1, refpos2, svpos1, svpos2, max_match
    return case, refpos1, refpos2, svpos1, svpos2, max_match

def process (key):
    global sv
    global ref
    sv = SV[key]
    ref = REF[key]
    len1 = len(sv)
    len2 = len(ref)
    i = 0
    j = 0
    Over = False
    while i < len1:
        while j < len2:
            # print("No. " + key + " checking: i = " + str(i) + " j = " + str(j))
            if sv[i] == ref[j]:
                i += 1
                j += 1
            else:
                # find max match
                case, refpos1, refpos2, svpos1, svpos2, max_match = greedy_search(i, j, len1, len2)
                if case == 1:
                    # insert
                    fout.write ("INS " + key + " " + str(refpos1) + " " + str(refpos2) + "\n")
                    i = svpos2 + max_match
                    j = refpos1 + max_match
                    if i >= len1 or j >= len2:
                        Over = True
                        break
                elif case == 2:
                    # duplicate
                    fout.write ("DUP " + key + " " + str(2*refpos1 - refpos2) + " " + str(refpos1) + "\n")
                    i = svpos2 + max_match
                    j = refpos1 + max_match
                    if i >= len1 or j >= len2:
                        Over = True
                        break
                elif case == 3:
                    # delete
                    fout.write ("DEL " + key + " " + str(refpos1) + " " + str(refpos2) + "\n")
                    i = svpos1 + max_match
                    j = refpos2 + max_match
                    if i >= len1 or j >= len2:
                        Over = True
                        break
                elif case == 4:
                    # inverse
                    fout.write ("INV " + key + " " + str(refpos1) + " " + str(refpos2) + "\n")
                    i = svpos2 + max_match
                    j = refpos2 + max_match
                    if i >= len1 or j >= len2:
                        print("over")
                        Over = True
                        break
                elif case == 5:
                    # translocation
                    tra = {}
                    tra['key'] = key
                    tra['refpos1'] = refpos1
                    tra['refpos2'] = refpos2
                    tra['svpos1'] = svpos1
                    tra['svpos2'] = svpos2
                    TRA.append(tra)
                    i = svpos2 + max_match
                    j = refpos2 + max_match
                if case == 0:
                    i += 1
                    j += 1
        if Over == True:
            break

start_time = time.time()
for key in SV.keys():
    process(key)

for tra1 in TRA:
    for tra2 in TRA:
        if tra2 == tra1 or tra1['key'] == '-1' or tra2['key'] == '-1':
            continue
        elif hash(REF[tra1['key']][tra1['refpos1']: tra1['refpos2']]) == hash(SV[tra2['key']][tra2['svpos1']: tra2['svpos2']]) and hash(REF[tra2['key']][tra2['refpos1']: tra2['refpos2']]) == hash(SV[tra1['key']][tra1['svpos1']: tra1['svpos2']]):
            fout.write("TRA " + tra1['key'] + " " + str(tra1['refpos1']) + " " + str(tra1['refpos2']) + tra2['key'] + " " + str(tra2['refpos1']) + " " + str(tra2['refpos2']) + "\n")
            tra2['key'] = '-1'
    tra1['key'] = '-1'
end_time = time.time()
print("Running time: " + str(end_time - start_time) + " s")
fout.close