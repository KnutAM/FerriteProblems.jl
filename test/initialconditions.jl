
@testset "Initial conditions" begin
    # Helper functions
    typebase(T::DataType) = T.name.wrapper
    typebase(v) = typebase(typeof(v))
    change_ip_order(ip::Interpolation, _::Nothing) = ip
    function change_ip_order(ip::Interpolation, order::Int)
        B = typebase(ip)
        Dim = Ferrite.getdim(ip)
        RefShape = Ferrite.getrefshape(ip)
        return B{Dim,RefShape,order}()
    end
    function test_case(;udim=1, CellType=Quadrilateral, iporder=nothing)
        grid = generate_grid(CellType, (10,10));
        ip_default = Ferrite.default_interpolation(CellType)
        ip = change_ip_order(ip_default, iporder)
        dh = DofHandler(grid); push!(dh, :u, udim, ip); push!(dh, :p, 1, ip); close!(dh);
        return dh
    end

    # Test cases
    for udim in (1,2)
        f = udim==1 ? Returns(1.0) : Returns(ones(Vec{udim}))
        for CellType in (Quadrilateral, Triangle, QuadraticQuadrilateral, QuadraticTriangle)
            for iporder in (1, 2)
                @testset "udim=$udim, cell=$CellType, ip order=$iporder" begin    
                    dh = test_case(;udim=udim, CellType=CellType, iporder=iporder)
                    a = zeros(ndofs(dh))
                    initial_conditions!(a, dh, :u, f)
                    @test sum(a)/length(a) ≈ udim/(1+udim)
                end
            end
        end
    end
end