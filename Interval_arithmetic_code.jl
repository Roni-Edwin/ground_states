using IntervalArithmetic


function phi_h(x)
    (pi..pi)*exp(-sqrt(2..2)*(pi..pi)*abs(x))*cos(sqrt(2..2)*(pi..pi)*abs(x)-(pi..pi)/4)
end

function T(x,n)
    (((4..4)*x^2*(n^2+(2..2)*x^2))/((n^2-(2..2)*x^2)^2))*(phi_h(x)-phi_h(n/sqrt(2..2)))
end

function deriv(k,n)
    (pi..pi)*(-2*(pi..pi))^k*cos(((pi..pi)*(k+1))/4)*(-1)^n*exp(-(pi..pi)*n)
end

function fanc_T(x,n)
    val=0
    for k=2:9
        val += (deriv(k,n)/factorial(k))*((x-n/sqrt(2..2))^(k-2))
    end
    val*((2..2)*x^2*(n^2+(2..2)*x^2))/((n+sqrt(2..2)*x)^2)
end

function I(x,n)
    if ((abs(x-n/sqrt(2..2)) ∩ (1..1)/100)==∅) & (abs(x-n/sqrt(2..2))<(1..1)/100)
        fanc_T(x,n)
    else
        T(x,n)
    end
end

function tilde_R(x)
    inter=0.0
    for n=1:50000
        inter += I(x,n)
    end
    inter+phi_h(x)-phi_h(0)
end

function tilde_R2(x)
    inter=0.0
    for n=2:50000
        inter += I(x,n)
    end
    inter+ (exp(-(pi..pi))*(pi..pi)^3)/sqrt(2..2) +phi_h(x)-phi_h(0)
end


pi_n=pi..pi
q2=sqrt(2..2)


lim=-(pi_n/q2)*tanh(pi_n/2)+(681..681)/(2^20)+(1..1)/(10^(10))
lim2=-(pi_n/q2)*tanh(pi_n/2)+(681..681)/(2^20)+(1..1)/(10^(10))+(1..1)/5

function check(lower,upper)
    beg=lower
    print(beg, "\n")
    last=16.0
    temp=beg..last
    while beg<upper+1/1024
        res=tilde_R(temp)
        iters=0
        while !(((res ∩ lim)==∅) & (res>lim))
            last=(beg+last)/2
            temp=beg..last
            res=tilde_R(temp)
            iters=iters+1
            if iters%10==0
                print("On iteration ",iters, "\n")
            end
            if iters>=20
                last=nextfloat(beg)
                iters=0
            end
        end
        beg=last
        last=beg+1
        temp=beg..last
        print("Checked up to ",beg, "\n")  
    end
    print("Checked on interval ", "[",lower,",",upper,"]","\n")
end

check(0,(693..693)/1000)
pats=[69,70,71,72]
for i in pats
    inp=((i..i)/100)..(((i+1)..(i+1))/100)
    res=tilde_R2(inp)
    if (((res ∩ lim2)==∅) & (res>lim2))
        print("Checked from 0.69 up to ",(i+1)/100, "\n")
    end
end
check((72..72)/100,1)
check(1,14)
print("\n", "In conclusion, inequality holds on interval [0,14]")

