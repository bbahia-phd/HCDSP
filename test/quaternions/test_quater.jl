# These tests were taken from the official Quaternions.jl package.
println("Testing the Quaternion type")
@test isa(typeof(Quaternion),DataType)
@test Quaternion <: Number
println("... tests passed!")

println("Testing Quaternion Constructor")
a = Quaternion(1)
b = Quaternion(1+im)
c = Quaternion(1+im,2(1+im))
q = Quaternion(1,1,1)

@test isa(a,Quaternion)
@test isa(b,Quaternion)
@test isa(c,Quaternion)
@test isa(q,Quaternion)
println("... tests passed!")

println("Testing quaternion alias")
a = quaternion(1)
b = quaternion(1+im)
c = quaternion(1+im,2(1+im))
q = quaternion(1,1,1)

@test isa(a,Quaternion)
@test isa(b,Quaternion)
@test isa(c,Quaternion)
@test isa(q,Quaternion)
println("... tests passed!")

println("Testing quater alias")
a = quater(1)
b = quater(1+im)
c = quater(1+im,2(1+im))
q = quater(1,1,1)

@test isa(a,Quaternion)
@test isa(b,Quaternion)
@test isa(c,Quaternion)
@test isa(q,Quaternion)
println("... tests passed!")


println("Testing type conversion")
a = convert(Quaternion,randn())
b = convert(Quaternion, randn() + randn()*im)
c = convert(Quaternion{Float64},randn())
q = convert(Float64,a)

@test isa(a,Quaternion)
@test isa(b,Quaternion)
@test isa(c,Quaternion)
@test isa(q,Float64)
println("... tests passed!")

println("Testing constructor on arrays")
A = zeros(10,10); B = zeros(10,10); C = zeros(10,10); D = zeros(10,10)

Q1 = quaternion(complex(A),complex(C))
@test isa(Q1,Array{Quaternion{Float64},2})

Q2 = quaternion(A,B,C,D)
@test isa(Q2,Array{Quaternion{Float64},2})

@test Q1 == Q2

Q = quaternion(A,B,C)
@test isa(Q,Array{Quaternion{Float64},2})

A = quaternion(A)
@test isa(A,Array{Quaternion{Float64},2})

B = quaternion(B + im .* C)
@test isa(B,Array{Quaternion{Float64},2})

println("... tests passed!")

println("Testing quaternion algebra")
test_associative(x, y, z, ⊗) = @test (x ⊗ y) ⊗ z ≈ x ⊗ (y ⊗ z)
test_commutative(x, y, ⊗) = @test x ⊗ y ≈ y ⊗ x
test_inverse(x, eins, ⊗, inv) = (@test x ⊗ inv(x) ≈ eins; @test inv(x) ⊗ x ≈ eins)
test_neutral(x, eins, ⊗) = (@test x ⊗ eins ≈ x; @test eins ⊗ x ≈ x)
test_monoid(x, y, z, ⊗, eins) = (test_associative(x, y, z, ⊗); test_neutral(x, eins, ⊗))
test_group(x, y, z, ⊗, eins, inv) = (test_monoid(x, y, z, ⊗, eins);test_inverse(x, eins, ⊗, inv))
test_multiplicative(x, y, ⊗, f) = @test f(x ⊗ y) ≈ f(x) ⊗ f(y)

# creating random examples
println("Defining random samples ... ")
sample(QT::Type{Quaternion{T}}) where {T <: Integer} = QT(rand(-100:100, 4)...)
sample(QT::Type{Quaternion{T}}) where {T <: AbstractFloat} = rand(Bool) ? randq() : nrandq()

sample(CT::Type{Complex{T}}) where {T <: Integer} = CT(rand(-100:100, 2)...)
sample(CT::Type{Complex{T}}) where {T <: AbstractFloat} = CT(randn(2)...)

sample(T, n) = T[sample(T) for _ in 1:n]

println("Testing algebraic properties of quaternions ... ")
# test algebraic properties of quaternions
for _ in 1:10, T in (Float32, Float64, Int32, Int64)
    q, q1, q2, q3 = sample(Quaternion{T}, 4)

    # skewfield
    test_group(q1, q2, q3, +, zero(q), -)
    test_group(q1, q2, q3, *, one(q), inv)
    test_multiplicative(q1, q2, *, norm)

    # complex embedding
    c1, c2 = sample(Complex{T}, 2)
    test_multiplicative(c1, c2, *, Quaternion)
    test_multiplicative(c1, c2, +, Quaternion)
end
println("... tests passed!")

println("Testing quaternion normalization ... ")
# test normalize
let
   q = randq()
   @test norm(normalize(q)) ≈ 1
   @test !isunit(q)
   @test q ≈ norm(q) * normalize(q)
   qn = nrandq()
   @test isapprox(norm(qn),1)
   @test isapprox(normalize(qn),qn) 
end
println("... tests passed!")


println("Testing quaternion involution ... ")
# test involutions
let
   q = randq()
   @test invi(q) == involution(q,quaternion(0,1,0,0))
   @test invj(q) == involution(q,quaternion(0,0,1,0))
   @test invk(q) == involution(q,quaternion(0,0,0,1))
end
println("... tests passed!")

println("Testing quaternion functions ... ")
# test overwritten math functions
for _ in 1:100
    let 
        c = Complex(randn(2)...)
        q, q2 = sample(Quaternion{Float64}, 4)
        unary_funs = [exp, log, sin, cos, sqrt, inv, conj, abs2, norm]

        for fun in unary_funs
            @test fun(Quaternion(c)) ≈ Quaternion(fun(c))
            @test q2 * fun(q) * inv(q2) ≈ fun(q2 * q * inv(q2))
        end

        @test exp(log(q)) ≈ q
        @test exp(zero(q)) ≈ one(q)
    end
end
println("... tests passed!")
