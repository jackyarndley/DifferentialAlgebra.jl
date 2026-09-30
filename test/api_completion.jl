using Test, LinearAlgebra, SpecialFunctions, TOML

@testset "Independent C-core coefficients" begin
    DifferentialAlgebra.init(6,3)
    x,y,z = DifferentialAlgebra.identity()
    p = 2+.1x+.2y+.05x*z+.03y^2
    data = TOML.parsefile(joinpath(@__DIR__,"data","c_core_coefficients.toml"))
    indices = DifferentialAlgebra.getMultiIndices(6,3)
    functions = ("SquareRoot"=>sqrt,"MultiplicativeInverse"=>inv,"Exponential"=>exp,
        "Logarithm"=>log,"Sine"=>sin,"Cosine"=>cos,"Tangent"=>tan,"ArcTangent"=>atan,
        "HyperbolicTangent"=>tanh,"GammaFunction"=>gamma,"LogGammaFunction"=>loggamma,
        "ErrorFunction"=>erf)
    for (name,f) in functions
        result = f(p)
        actual = [DifferentialAlgebra.getCoefficient(result,jj) for jj in indices]
        @test maximum(abs.(actual-data[name]["coefficients"])) < 1e-12
    end
end

# Every matrix access checks the public order, including accesses inside the
# eigenpair lifting loop. This deterministically detects hidden setTO calls.
struct OrderCheckedMatrix{T} <: AbstractMatrix{T}
    data::Matrix{T}
    order::Int
end
Base.size(A::OrderCheckedMatrix) = size(A.data)
function Base.getindex(A::OrderCheckedMatrix,i::Int,j::Int)
    DifferentialAlgebra.getTO() == A.order || error("eigh changed the global truncation order")
    A.data[i,j]
end

@testset "API audit regressions" begin
    for T in (Float32,Float64,BigFloat), budget in (0,32*1024^2)
        DifferentialAlgebra.init(6,3;table_bytes=budget)
        x,y,z = DifferentialAlgebra.identity(T)
        p = 3+x+2y+x^2*y+3x*z^2+x^5
        tolerance = T === Float32 ? T(2e-5) : T === Float64 ? T(2e-13) : T(2)^(-200)
        @test DifferentialAlgebra.norm(DifferentialAlgebra.plug(p,1,2)-(5+6y+6z^2+32)) == 0
        @test DifferentialAlgebra.norm(DifferentialAlgebra.plug(p,1)-(3+2y)) == 0
        @test DifferentialAlgebra.norm(DifferentialAlgebra.replaceVariable(p,1,2,2)-(3+4y+4y^3+6y*z^2+32y^5)) == 0
        for v in 1:3
            args = [x,y,z]
            args[v] = 2args[v]+1
            @test DifferentialAlgebra.norm(DifferentialAlgebra.translateVariable(p,v,2,1)-DifferentialAlgebra.evaluate(p,args)) < tolerance
            args[v] = -[x,y,z][v]
            @test DifferentialAlgebra.norm(DifferentialAlgebra.scaleVariable(p,v,-1)-DifferentialAlgebra.evaluate(p,args)) < tolerance
        end
        @test_throws ArgumentError DifferentialAlgebra.plug(p,0)
        @test_throws ArgumentError DifferentialAlgebra.replaceVariable(p,1,4)
        @test DifferentialAlgebra.evalMonomials(p,2x+3y+x^2*y) == T(9)
        @test DifferentialAlgebra.norm(DifferentialAlgebra.filterMonomials(p,x+y)-(x+2y)) == 0
        @test DifferentialAlgebra.nterms(p) == 6
        @test DifferentialAlgebra.maxNorm(p) == 3
        @test DifferentialAlgebra.norm(DifferentialAlgebra.monomial([2,1,0],T(3))-3x^2*y) == 0
        @test DifferentialAlgebra.nterms(DifferentialAlgebra.filled(T(2))) == DifferentialAlgebra.getMaxMonomials()
        @test DifferentialAlgebra.filled(T(2)) isa DA{T}
        q = inv(1-x/4)
        estimates,errors = DifferentialAlgebra.estimNorm(q,0,0,8;errors=true)
        @test maximum(abs.(estimates-[T(1)/4^i for i in 0:8])) < tolerance
        @test all(>=(0),errors) && length(errors) == 7
        @test DifferentialAlgebra.convRadius(q,T(1)/2^14) ≈ one(T) rtol=tolerance
        @test_throws ArgumentError DifferentialAlgebra.convRadius(q,-1)
        @test DifferentialAlgebra.isInitialized()
        DifferentialAlgebra.setTO(2); DifferentialAlgebra.pushTO()
        @test DifferentialAlgebra.getTO() == 6
        DifferentialAlgebra.popTO()
        @test DifferentialAlgebra.getTO() == 2
        DifferentialAlgebra.setTO()
        @test DifferentialAlgebra.getTO() == 6
        @test DifferentialAlgebra.norm((p/q)*q-p) < tolerance
        @test DifferentialAlgebra.norm(muladd(T(2),p,q)-(2p+q)) < tolerance
        @test DifferentialAlgebra.norm(muladd(p,T(2),q)-(2p+q)) < tolerance
        @test DifferentialAlgebra.norm(muladd(p,q,p)-(p*q+p)) < tolerance
        @test DifferentialAlgebra.norm(sqrt(2+x+y*z)^2-(2+x+y*z)) < tolerance
        @test DifferentialAlgebra.norm(DifferentialAlgebra.fromString(DifferentialAlgebra.toString(p),T)-p) == 0
        @test parse(DA{T},DifferentialAlgebra.toString(q)) isa DA{T}
        @test DifferentialAlgebra.norm(parse(DA{T},DifferentialAlgebra.toString(q))-q) == 0
        @test iszero(DifferentialAlgebra.fromString(DifferentialAlgebra.toString(zero(p)),T))
        @test asin(DA{T}(1)) isa DA{T}
        @test DifferentialAlgebra.cons(acos(DA{T}(-1))) == acos(-one(T))
        @test iszero(acosh(DA{T}(1)))
        @test_throws DifferentialAlgebra.DAError asin(1+x)
        v = AlgebraicVector([2+x,3+y])
        for f in (sqrt,exp,log,sinh,cosh,tanh,asinh,acosh,inv)
            @test maximum(DifferentialAlgebra.norm,f(v)-AlgebraicVector(f.(v.data))) < tolerance
        end
        @test maximum(DifferentialAlgebra.norm,v^2-DifferentialAlgebra.sqr(v)) == 0
        @test DifferentialAlgebra.norm(DifferentialAlgebra.vnorm(normalize(v))-1) < 10tolerance
        @test length(AlgebraicVector{DA{T}}()) == 0
        filled = AlgebraicVector{DA{T}}(2,x)
        DifferentialAlgebra.setCoefficient!(filled[1],[1,0,0],3)
        @test DifferentialAlgebra.getCoefficient(filled[2],[1,0,0]) == 1
        A = OrderCheckedMatrix([2+x/10 one(x)/5;one(x)/5 4+y/10],6)
        values,vectors = DifferentialAlgebra.eigh(A)
        @test maximum(DifferentialAlgebra.norm,A*vectors-vectors*Diagonal(values)) < 20tolerance
        @test maximum(DifferentialAlgebra.norm,vectors'*vectors-I) < 20tolerance
        @test isempty(first(DifferentialAlgebra.eigh(Matrix{DA{T}}(undef,0,0))))
        @test size(DifferentialAlgebra.linear(DA{T}[])) == (0,3)
        @test size(DifferentialAlgebra.hess_stack(DA{T}[])) == (3,3,0)
    end

    DifferentialAlgebra.init(5,2)
    x,y = DifferentialAlgebra.identity()
    for T in (Float32,Float64,BigFloat)
        v = AlgebraicVector(T[1,4,9])
        @test DifferentialAlgebra.sqr(v).data == T[1,16,81]
        @test DifferentialAlgebra.root(v).data ≈ T[1,2,3]
        @test DifferentialAlgebra.isrt(v).data ≈ T[1,one(T)/2,one(T)/3]
        @test DifferentialAlgebra.icrt(AlgebraicVector(T[-1,-8])).data ≈ T[-1,-1/2]
    end
    @test DifferentialAlgebra.cons(hypot(DA(3e200),DA(4e200))) ≈ 5e200
    @test DifferentialAlgebra.cons(hypot(DA(3e-200),DA(4e-200))) ≈ 5e-200
    @test DifferentialAlgebra.cons(round(DA(Inf))) == Inf
    @test DifferentialAlgebra.cons(trunc(DA(-Inf))) == -Inf
    @test DifferentialAlgebra.norm(mod(-1.5+x,2)-(0.5+x)) == 0
    setprecision(256) do
        DifferentialAlgebra.setEps(BigFloat(1)-BigFloat(2)^(-100))
        @test DifferentialAlgebra.cons(DA(1.0)) == 1
        @test DifferentialAlgebra.cons(DA(1f0)) == 1
        DifferentialAlgebra.setEps(0)
    end
    @test DifferentialAlgebra.getCoefficient((1e308x+1e-308y)^2,[1,1]) ≈ 2
    # Hand-specified C-core blob for 1.5 + 2*x*y at order five, two variables.
    fixture = hex2bytes("4441301e05000000010000000100000002000000" *
        "0000000000000000000000000000f83f01000000010000000000000000000040")
    @test DifferentialAlgebra.norm(read(IOBuffer(fixture),DA)-(1.5+2x*y)) == 0
    io = IOBuffer()
    @test write(io,1.5+2x*y) == length(fixture)
    @test take!(io) == fixture
    for p in (zero(x),x,1+x+2y+x*y+x^5)
        io = IOBuffer(); write(io,p); seekstart(io)
        @test DifferentialAlgebra.norm(read(io,DA)-p) == 0
    end
    @test_throws ArgumentError read(IOBuffer(zeros(UInt8,36)),DA)
    @test_throws EOFError read(IOBuffer(fixture[1:25]),DA)
    @test_throws ArgumentError write(IOBuffer(),DA{BigFloat}(1))
    @test_throws ArgumentError DifferentialAlgebra.fromString("incorrect header")
    @test_throws ArgumentError DifferentialAlgebra.fromString("I COEFFICIENT ORDER EXPONENTS\n1 2.0 3 1 1")
    text = "     I  COEFFICIENT              ORDER EXPONENTS\n" *
        "1 1.5D+00 0 0 0\n2 2.0D+00 2 1 1\n------------------------------------------------\n"
    @test DifferentialAlgebra.norm(DifferentialAlgebra.fromString(text)-(1.5+2x*y)) == 0
    DifferentialAlgebra.init(1,1)
    @test DifferentialAlgebra.cons(read(IOBuffer(fixture),DA)) == 1.5
    @test_throws ArgumentError DifferentialAlgebra.getMultiIndices(1000,1000)
end

@testset "Complete special functions at coefficient precision" begin
    for bits in (128,256)
        setprecision(bits) do
            DifferentialAlgebra.init(4,1)
            x = DifferentialAlgebra.variable(1,BigFloat)
            tol = BigFloat(2)^(-bits+30)
            for n in 1:5
                expected = (-1)^(n+1)*factorial(big(n))*(zeta(BigFloat(n+1))-1)
                @test DifferentialAlgebra.cons(polygamma(n,2+zero(x))) ≈ expected rtol=tol
                a = -BigFloat(1)/2+x/10
                @test DifferentialAlgebra.norm(polygamma(n,a+1)-polygamma(n,a)-(-1)^n*factorial(big(n))*a^(-n-1)) < 100tol
            end
            # Independent 100-decimal-digit mpmath 1.4.1 reference values.
            for (f,expected) in (
                (p->besseli(2,p),"0.809196566007602118254159830538643220568413967382672328252350891607223341316661635156423216783493928"),
                (p->besselk(2,p),"0.2096093039998351971497158542154000893596098104942351105661592886876975266704431100028063620394724142"))
                @test DifferentialAlgebra.cons(f(BigFloat(17)/8+zero(x))) ≈ parse(BigFloat,expected) rtol=tol
            end
            expected = "344.6756950613542964686513088383632754343010336756080598637557812072204830205585200518156277649550717"
            @test DifferentialAlgebra.cons(polygamma(3,-BigFloat(5)/8+zero(x))) ≈ parse(BigFloat,expected) rtol=tol
            for center in (BigFloat(2),-BigFloat(1)/2)
                p = center+x/10
                @test gamma(p) isa DA{BigFloat}
                @test DifferentialAlgebra.norm(gamma(p+1)-p*gamma(p)) < 10tol
                @test abs(DifferentialAlgebra.evaluate(gamma(p),[BigFloat(1)/1000])-gamma(center+BigFloat(1)/10000)) < BigFloat(1e-16)
            end
            @test_throws DomainError polygamma(2,zero(x))
            for center in (BigFloat(1)/8,BigFloat(2),BigFloat(20),BigFloat(300)), n in (0,2)
                p = center+x/100
                i,k = besselix(n,p),besselkx(n,p)
                @test i isa DA{BigFloat} && k isa DA{BigFloat}
                @test DifferentialAlgebra.norm(i*besselkx(n+1,p)+besselix(n+1,p)*k-inv(p)) < 100tol
                @test DifferentialAlgebra.norm(besseli(n,p)*exp(-p)-i) < 100tol
                @test DifferentialAlgebra.norm(besselk(n,p)*exp(p)-k) < 100tol
            end
            @test DifferentialAlgebra.norm(besseli(1,-(2+x))+besseli(1,2+x)) < tol
            @test DifferentialAlgebra.norm(besseli(0,x)-(1+x^2/4+x^4/64)) < tol
            @test_throws DifferentialAlgebra.DAError besselix(0,x)
            @test_throws DomainError besselk(0,-one(x))
        end
    end
    DifferentialAlgebra.init(4,1)
    x = DifferentialAlgebra.variable(1,Float32)
    for f in (gamma,loggamma,digamma,p->polygamma(2,p),p->besseli(2,p),p->besselk(2,p),p->besselix(2,p),p->besselkx(2,p))
        @test f(2+x/10) isa DA{Float32}
    end
end
