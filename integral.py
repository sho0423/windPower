import math
import random
import numpy as np
import matplotlib.pyplot as plt

# ワイブル分布

class Weibull():
    def __init__(self):
        pass
        
    # 密度関数f, a:ワイブル係数, b: 尺度パラメータ
    def f(self, x, a, b):
        return (a/b) * ( (x/b)**(a - 1) ) * math.exp(-1 * (x/b)**a )
    
    # 分布関数F
    def F(self, x, a, b):
        return 1 - math.exp(-1 * (x/b)**a )
    
    # Fの逆関数
    def invF(self, y, a, b):
        return b * ( -1 * math.log(1 - y) )**(1/a)
    
    # サンプリング
    def sampling(self, num, a, b):
        sample=[]
        for i in range(num):
            rand = random.random()
            sample.append(self.invF(rand, a, b))
        return np.array(sample)
    
    # 散布図
    def scatter_pdf(self, ax, data, a, b, alpha=0.4):
        pdf=[]
        for d in data:
            pdf.append(self.f(d, a, b))

        ax.scatter(data, pdf, alpha=alpha)
        ax.set_title("Weibull (alpha={:.2f}, beta={:.2f})".format(a, b))
        ax.set_xlabel("value")
        ax.set_ylabel("probability density")
        ax.set_ylim(0, 1)
        
    # グラフ
    def plot_pdf(self, ax, data, a, b):
        pdf=[]
        for d in data:
            pdf.append(self.f(d, a, b))

        ax.plot(data, pdf)
        ax.set_title("Weibull (alpha={:.2f}, beta={:.2f})".format(a, b))
        ax.set_xlabel("value")
        ax.set_ylabel("probability density")
        ax.set_ylim(0, 0.25)

data = np.linspace(0, 20, 100)
'''ワイブル分布プロット用
fig, ax = plt.subplots(3, 3, figsize=(18, 18))

Weibull().plot_pdf(ax[0, 0], data, 1, 1)
Weibull().plot_pdf(ax[0, 1], data, 1, 2)
Weibull().plot_pdf(ax[0, 2], data, 1, 3)
Weibull().plot_pdf(ax[1, 0], data, 2, 1)
Weibull().plot_pdf(ax[1, 1], data, 2, math.sqrt(2)*5.585)
Weibull().plot_pdf(ax[1, 2], data, 2, 3)
Weibull().plot_pdf(ax[2, 0], data, 3, 1)
Weibull().plot_pdf(ax[2, 1], data, 3, 2)
Weibull().plot_pdf(ax[2, 2], data, 3, 3)

plt.savefig("out/Weibull_Sample.png")
'''

sigma = 7*math.sqrt(2/math.pi)

def v_prob(v):
  return Weibull().f(v, 2, math.sqrt(2)*sigma)

def p(v):
  if v<=3 or v>= 25:
    return 0
  elif v>= 11 and v <25:
    return 2000
  else:
    return 2000*((v-3)/8)**3

def prod(v):  #被積分函数
  return p(v)*v_prob(v)

def integ(func,acc,st,end):  #積分
  s = 0
  v = st
  d = (end - st)/acc
  for i in range(acc):
    s += func(v)*d
    v += d

  return s

print(integ(prod,10**6,0,25))