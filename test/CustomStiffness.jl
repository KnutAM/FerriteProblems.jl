import FerriteProblems: CustomStiffness, get_stiffness_type
# Have this wrapper to not directly modify the behavior for HeatEquation
struct CSTestWrap{M}
    m::M
end
function FerriteAssembly.element_routine!(Ke, re, state, ae, m::CSTestWrap, args...; kwargs...)
    FerriteAssembly.element_routine!(Ke, re, state, ae, m.m, args...; kwargs...)
end
function FerriteAssembly.element_residual!(re, state, ae, m::CSTestWrap, args...; kwargs...)
    FerriteAssembly.element_residual!(re, state, ae, m.m, args...; kwargs...)
end

function FerriteAssembly.element_routine!(Ke, re, state, ae, m::CustomStiffness{<:CSTestWrap}, args...; kwargs...)
    FerriteAssembly.element_routine!(Ke, re, state, ae, m.material, args...; kwargs...)
    if get_stiffness_type(m) == :negative
        for i in 1:size(Ke,1)
            Ke[i,i] *= -1
        end
    end
end

@testset "CustomStiffness" begin
    grid = generate_grid(Quadrilateral, (2,2))
    dh = DofHandler(grid)
    ip = Lagrange{RefQuadrilateral,1}()
    add!(dh, :u, ip)
    close!(dh)
    cv = CellValues(QuadratureRule{RefQuadrilateral}(2), ip)
    mw = CSTestWrap(FerriteAssembly.ExampleElements.StationaryFourier(1.0))
    m = CustomStiffness(mw, :positive)
    buffer = setup_domainbuffer(DomainSpec(dh, m, cv))
    buffer_ad = setup_domainbuffer(DomainSpec(dh, m, cv); autodiffbuffer=true)
    K = create_sparsity_pattern(dh)
    r = zeros(ndofs(dh))
    a = rand(ndofs(dh))

    assembler = start_assemble(K, r)
    work!(assembler, buffer)
    d1 = sum(diag(collect(K)))
    m.stiffness_type = :negative
    assembler = start_assemble(K, r)
    work!(assembler, buffer)
    d2 = sum(diag(collect(K)))
    @test d1 ≈ -d2

    assembler = start_assemble(K, r)
    work!(assembler, buffer_ad; a=a)
    @test d2 ≈ sum(diag(collect(K)))

    r2 = rand(length(r))
    r_assembler = ReAssembler(r2)
    work!(r_assembler, buffer; a=a)
    @test r ≈ r2

    # Check that correct error is thrown if not implemented
    m_throw = CustomStiffness(mw.m, :positive) # element_routine! not implemented 
    b_throw = setup_domainbuffer(DomainSpec(dh, m_throw, cv))
    assembler = start_assemble(K, r)
    @test_throws ArgumentError work!(assembler, b_throw)
    msg = "You must implement element_routine!, normally for m::CustomStiffness{<:"
    @test_throws msg work!(assembler, b_throw)
end