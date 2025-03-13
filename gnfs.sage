from math import isqrt, log, sqrt, exp, floor, ceil

def smoothTest(x,factorBase):
    xAbs=abs(x)
    exponents=[0]*len(factorBase)
    for i,p in enumerate(factorBase):
        while xAbs % p==0:
            exponents[i]+=1
            xAbs//=p
        if xAbs==1:
            break
    return exponents if xAbs==1 else None
def algorithm(n):
    if n<2:
        return n
    if is_prime(n):
        return n
    d=3
    m=Integer(floor(n**(1/d)))
    R=PolynomialRing(ZZ,'x')
    x=R.gen()
    poly=x**d-m**d+n
    K=NumberField(poly,'a')
    a=K.gen()
    B=ceil(exp(sqrt(log(n)*log(log(n)))))
    factorBase=list(prime_range(2,B))
    sieveBound=B
    logFactorBase=[log(p) for p in factorBase]
    relations=[]
    relationData=[]
    threshold=len(factorBase)*2+5
    for bVal in range(1,sieveBound+1):
        candidates={}
        for aVal in range(-sieveBound,sieveBound+1):
            A=aVal-bVal*m
            if A==0:
                continue
            candidates[aVal]=log(abs(A))
        for i,p in enumerate(factorBase):
            r=(bVal*m) % p
            for aVal in range(-sieveBound,sieveBound+1):
                if (aVal-r) % p==0 and aVal in candidates:
                    candidates[aVal]-=logFactorBase[i]
        for aVal,rem in candidates.items():
            if rem<1:
                ratVal=aVal-bVal*m
                sRat=smoothTest(ratVal,factorBase)
                if sRat is None:
                    continue
                elem=aVal-bVal*a
                algVal=elem.norm()
                sAlg=smoothTest(algVal,factorBase)
                if sAlg is None:
                    continue
                relations.append([e%2 for e in sRat]+[e%2 for e in sAlg])
                relationData.append((ratVal,algVal))
        if len(relations)>=threshold:
            M=Matrix(GF(2),relations)
            ker=M.right_kernel()
            if ker.dimension()>0:
                dep=ker.basis()[0]
                xProd=1
                yProd=1
                for i,coeff in enumerate(dep):
                    if coeff:
                        xProd*=relationData[i][0]
                        yProd*=relationData[i][1]
                if xProd<0:
                    xProd=-xProd
                if yProd<0:
                    yProd=-yProd
                xRoot=isqrt(xProd)
                yRoot=isqrt(yProd)
                factorCandidate=gcd(xRoot-yRoot,n)
                if factorCandidate not in [1,n]:
                    return factorCandidate
    return None

print(algorithm(3*5))
