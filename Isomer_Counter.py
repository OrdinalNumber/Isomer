import numpy as np
from rdkit import Chem
import re
import time

s = input("分子式を入力してください : ")

start = time.time()

def count_elem(exp):
  d = {}
  m = re.findall("[A-Z][a-z]*[0-9]*|\(.+?\)[0-9]*",exp)
  for i in m:
    if i[0] == "(":
      t = re.match("(\(.+?\))([0-9]*)",i).groups()
      d_ = count_elem(t[0][1:-1])
      for k in d_:
        if k not in d:
          d[k] = 0
        if t[1] == "":
          d[k] += d_[k]
        else:
          d[k] += d_[k]*int(t[1])
    else:
      t = re.match("([A-Z][a-z]*)([0-9]*)",i).groups()
      if t[1] == "":
        if t[0] not in d:
          d[t[0]] = 0
        d[t[0]] += 1
      else:
        if t[0] not in d:
          d[t[0]] = 0
        d[t[0]] += int(t[1])
  return d

mols_dict = count_elem(s)
print(mols_dict)
if "H" not in mols_dict:
  mols_dict["H"] = 0


mols_list = sum([[k]*v for k,v in mols_dict.items() if k != "H"],[])
print(mols_list)

element_data = {
    'H': 1,
    'He': 0,
    'Li': 1,
    'Be': 2,
    'B': 3,
    'C': 4,
    'N': 3,
    'O': 2,
    'F': 1,
    'Ne': 0,
    'Na': 1,
    'Mg': 2,
    'Al': 3,
    'Si': 4,
    'P': 3,
    'S': 2,
    'Cl': 1,
    'Ar': 0,
    'K': 1,
    'Ca': 2,
    'Rb': 1,   # アルカリ金属: ルビジウム (Rubidium)
    'Sr': 2,   # アルカリ土類金属: ストロンチウム (Strontium)
    'Cs': 1,   # アルカリ金属: セシウム (Cesium)
    'Ba': 2,   # アルカリ土類金属: バリウム (Barium)
    'Fr': 1,   # アルカリ金属: フランシウム (Francium)
    'Ra': 2,   # アルカリ土類金属: ラジウム (Radium)
    'Cl': 1,   # ハロゲン: 塩素 (Chlorine)
    'Br': 1,   # ハロゲン: 臭素 (Bromine)
    'I': 1,    # ハロゲン: ヨウ素 (Iodine)
    'At': 1    # ハロゲン: アスタチン (Astatine)
}

mols_count = len(mols_list) #原子の個数
mols_upper = [element_data[i] for i in mols_list] #各原子の腕の本数
Iu = 1-mols_count-mols_dict["H"]+(sum(mols_upper)+mols_dict["H"])//2 #水素とくっつかない結合の本数*2
connection_count = sum(mols_upper)-mols_dict["H"]
filter = mols_list #定義：同じ数字or文字が連続しているカタマリで構成されている、同じ数字の塊が2度登場することはない。原子番号を表す。
print("==============")
print(connection_count)
print(mols_count)
print(mols_upper)
print(filter)
print("==============")
print(Iu)
print("==============")

'''
connection_count = 18 #水素とくっつかない結合の本数*2
mols_count = 6 #原子の個数
mols_upper = [4,4,4,4,4,4] #各原子の腕の本数
filter = [6,6,6,6,6,6] #定義：同じ数字が連続しているカタマリで構成されている、同じ数字の塊が2度登場することはない。原子番号を表す。
print("==================")
print(connection_count)
print(mols_count)
print(mols_upper)
print(filter)
print("==================")
'''

def div_hydro(connect,limit): #水素の分配(逆説的にね(水素以外との結合数を計算して))
  if len(limit) == 0 or connect <= 0:
    return []
  elif len(limit) == 1:
    if connect <= limit[0]:
      return [[connect]]
    else:
      return []
  if sum(limit) < connect:
    return []
  res = []
  for i in range(1,limit[0]+1):
    array = div_hydro(connect-i,limit[1:])
    for k in array:
      res.append([i]+k)
  return res

#print(div_hydro(connection_count,mols_upper))
X = div_hydro(connection_count,mols_upper)
filter_copy = filter
filter_array = []
c = 1
b = filter_copy[0]
while filter_copy:
  if filter_copy[c] == b:
    c += 1
  else:
    b = filter_copy[c]
    filter_array.append(c)
    filter_copy = filter_copy[c:]
    c = 1
  if c == len(filter_copy):
    filter_array.append(c)
    filter_copy = filter_copy[c:]
print(filter_array)
Y = []
for x in X:
  x_copy = x
  divided_x = []
  for i in filter_array:
    divided_x.append(x_copy[0:i])
    #print("in function",x_copy[0:i])
    x_copy = x_copy[i:]
  divided_x = [sorted(e) for e in divided_x]
  canonical_x = sum(divided_x,[])
  if not canonical_x in Y:
    Y.append(canonical_x)
Y.sort()
#print(len(Y),Y)


#完成版、ラプラシアン行列 to SMILES記法
def Laplacian2SMILES(sample,mols_arg):
  mols = list(mols_arg)
  point = 0
  num = 1
  stuck = [point]
  deepth = [0]
  path = [point]
  seen = [False]*len(mols)
  res_id = []
  res_list = []
  while stuck:
    point = stuck.pop()
    d = deepth.pop()
    res_id.append(point)
    path = path[:d]+[point]
    res_list.append(path)
    seen[point] = True
    for i in range(len(mols)):
      if sample[point][i] < 0:
        if not seen[i]:
          if i in stuck:
            index = stuck.index(i)
            del stuck[index]
            del deepth[index]
          stuck.append(i)
          deepth.append(d+1)
        elif i in path and i != path[-2]:
          mols[point] += f"{num}"
          if sample[point][i] ==-1:
            mols[i] += f"{num}"
          if sample[point][i] ==-2:
            mols[i] += f"={num}"
          if sample[point][i] ==-3:
            mols[i] += f"#{num}"
          num += 1
          
  #関数内部の関数定義
  def list2smiles(id_list,path_list,elem_list,deep = 0,before = -1):
    if len(id_list) != len(path_list):
      raise ValueError("do not much length 'id_list' and 'path_list'")
    if len(id_list) == 0:
      return ""
    result = ""
    M = max(len(path) for path in path_list)
    for i in range(M):
      array = []
      for k in range(len(id_list)):
        if len(path_list[k]) > i:
          array.append(path_list[k][i])
      connect_set = [e for e in array if e >= 0]
      connect_set = sorted(set(connect_set), key=connect_set.index)
      if len(connect_set) == 1:
        if deep == 0 and before == -1:
          result += elem_list[id_list[i]]
        else:
          if sample[before][id_list[i]] == -1:
            result += elem_list[id_list[i]]
          elif sample[before][id_list[i]] == -2:
            result += "=" + elem_list[id_list[i]]
          elif sample[before][id_list[i]] == -3:
            result += "#" + elem_list[id_list[i]]
        before = id_list[i]
      else:
        m = connect_set[-1]
        for n in connect_set:
          next_id_list = []
          next_path_list = []
          for k in range(len(id_list)):
            if len(path_list[k]) > i and path_list[k][i] == n:
              next_id_list.append(id_list[k])
              next_path_list.append(path_list[k][i:])
          if n != m:
            result += "(" + list2smiles(next_id_list,next_path_list,elem_list,before = before,deep=deep+1) + ")"
          else:
            result += list2smiles(next_id_list,next_path_list,elem_list,before = before,deep=deep+1)
        return result
    return result

  res = list2smiles(res_id,res_list,mols)
  return res

#======================================================================
partition_memo = {}

def partition(n, k):
  try:
    return partition_memo[(n,k)]
  except:
    array = partition_(n, k)
    partition_memo[(n,k)] = tuple(array)
    return partition_memo[(n,k)]

def partition_(n, k):
  if n < 0:
    return []
  if k == 1:
    return [[n]]
  res = []
  for i in range(n+1):
    res.extend([array + [i] for array in partition_(n - i, k - 1)])
  return res
#======================================================================

def generate_Lap(L,already_conection = None,deep = 0):
  if len(L) == 0:
    return []
  elif len(L) == 1:
    return [[[0]]]
  if already_conection == None:
    already_conection = [0]*len(L)
  result = []
  for x in partition(L[0]-already_conection[0],len(L)-1):
    for m in generate_Lap(L[1:],[a+b for a,b in zip(already_conection[1:],x)],deep=deep+1):
      res = [[0]+x]+[[x[i]]+l for i,l in enumerate(m)]
      result.append(res)
  if deep == 0:
    result_ = []
    A = np.zeros((len(L),len(L)))
    for i in range(len(L)):
      A[i][i] = L[i]
    for res in result:
      if int(sum(res[-1])) == int(L[-1]):
        #print(sum(res[-1]),L[-1],sum(res[-1]) == L[-1],len(result_))
        array = np.array(res)*(-1)+A
        #print(Laplacian2SMILES(array,mols_list))
        if np.linalg.matrix_rank(array) == mols_count-1:
          smile = Chem.MolToSmiles(Chem.MolFromSmiles(Laplacian2SMILES(array,mols_list)),isomericSmiles = False,kekuleSmiles = True)
          #time.sleep(0.2)
          #print(smile)
          #print(array)
          if smile not in result_:
            print(smile)
            result_.append(smile)
      #print(sum(res[-1]),L[-1],type(sum(res[-1])),type(L[-1]),sum(res[-1])==L[-1],int(sum(res[-1]))==int(L[-1]))
    return result_
  else:
    return result

K = []
for L in Y:
  Z = generate_Lap(L)
  print("partial",len(Z),L)
  K += Z
print("total",len(K))

#K = [array for array in K if np.linalg.matrix_rank(array) == mols_count-1]


R_ = K
R_ = sorted(set(R_), key=R_.index)
#R_ = [Chem.MolToSmiles(Chem.MolFromSmiles(x),isomericSmiles = False,kekuleSmiles = True) for x in R_]
#R_ = sorted(set(R_), key=R_.index)

end = time.time()

print("===================")
print("Isomer of ",s)
print("the number is :",len(R_))
print("processing time : ",f'{end-start}s')
print("SMILES:",R_)
print("===================")


res = f'Isomer of {s}\nthe number is : {len(R_)}\nprocessing time : {end-start}s\nSMILES:{R_}'
with open(f"Isomer_of_{s}.txt", mode='w') as f:
    f.write(res)

print("書き込み終了(「end」と入力すれば終了します)")

text = ""
while text != "end":
  text = input()