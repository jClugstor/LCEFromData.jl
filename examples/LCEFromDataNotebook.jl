### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# ╔═╡ c6b1ad98-60b8-11ef-12c9-7f6a2729e4b3
using Pkg 

# ╔═╡ 77862b58-3162-420c-b34e-183d9337b564
Pkg.add(url="https://github.com/jClugstor/LCEFromData.jl")

# ╔═╡ b44b2a30-d948-44b9-8804-8517eca4ce3c
using LCEFromData

# ╔═╡ 972d3b67-9b8d-46b0-83e1-540b76175516
using DynamicalSystems

# ╔═╡ 735bb737-de5d-44d5-aaa5-f31246b98b09
md"""
Because my package isn't in the general registry yet, you can just add it directly from my GitHub.
"""

# ╔═╡ 0a6ea318-7414-4943-95ba-b8144984d62d
md"# Example Data Generation"

# ╔═╡ 2cf1e8fe-7797-4dc8-85b4-c3a41f1445d6
md"""
We'll use DynamicalSystems.jl to generate some data to use.
"""

# ╔═╡ b874c254-222a-4f74-b924-cc785482189b
md"""
Just a function that represents the RHS of an ODE, e.g. $\dot{x} = f(x,t)$, this is $f(x,t)$. By convention Julia functions that mutate input objects have "!" appended to them. 
"""

# ╔═╡ 4cd5e33e-964c-46cb-96d5-a868a328bfcf
function LorenzSystem!(du, u, p, t)
        σ, β, ρ = p
        x, y, z = u
        du[1] = σ * (y - x)
        du[2] = x * (ρ - z) - y
        du[3] = x * y - β * z
end

# ╔═╡ bd78472f-de07-43ad-a21e-85ec6b8a6cce
md"""
This will set up the parameters and $u_0$. The `begin` block is just because this notebook environment is "reactive", so every cell should be one declaration, unless it's a `let` or `begin` block. For more details see [Pluto.jl](https://plutojl.org).
"""

# ╔═╡ ee2556f0-3c5e-4f85-92ac-37f7422c7057
begin
	pars = [16.0, 4.0, 45.92]
	u0 = [0.1, 0.1, 0.1]
end

# ╔═╡ 063da746-112c-4c6f-8427-606406870484
md"""
This uses sets up a DynamicalSystems.jl object to use `LorenzSystem!` as the "update rule" and `pars` and `u0` as parameters and initial condition. 
"""

# ╔═╡ ee89f54a-823c-4a78-9682-df4cad4349a2
lorenz_sys = CoupledODEs(LorenzSystem!,u0,pars)

# ╔═╡ a6f1252a-1518-46db-b52d-6c766f913186
Y,t = trajectory(lorenz_sys, 500.0, Δt = 0.1)

# ╔═╡ 23c26188-aa54-467b-99a6-be7f40a9655e
md"Then this gets the x state of the trajectory."

# ╔═╡ 7a26b34f-bb2d-4834-961b-52aadeb7b0d1
x_dat = Y[:,1] 

# ╔═╡ d24213c9-9328-421d-b3a1-cc00320c8d80
md"# LCEFromData"

# ╔═╡ 34571dbf-fc63-4a50-8ea5-9ff12a12e1ac
md"""First set up an `LCEProblem(timeseries, dt, dim, tau)`"""

# ╔═╡ ea6d6865-2470-481b-8d1b-4cd0e02468eb
LCE_prob = LCEProblem(x_dat, 0.1, 3, 3)

# ╔═╡ 306ad617-7e12-4582-a250-f4817004cc5c
md"Then construct an algorithm object with the desired parameters. In this case I'll just put in the `ks` parameter, the rest of the parameters have defaults (for all the keyword argument options you can see the docstring of `DivergenceAlgorithm` in the source code, or you can go to the Live Docs button on the bottom right and search for `DivergenceAlgorithm`."

# ╔═╡ 39ad81dd-391f-46e8-84b6-bfa8456d30b8
div_alg = DivergenceAlgorithm(0:1:100)

# ╔═╡ 2a7cee5a-9caa-4760-8923-e9362329f65a
md"Then you can solve the problem with the chosen algorithm using the `solve` function."

# ╔═╡ 80d53fb6-5052-473c-8f17-8ec4a5c88889
solve(LCE_prob, div_alg)

# ╔═╡ Cell order:
# ╠═c6b1ad98-60b8-11ef-12c9-7f6a2729e4b3
# ╟─735bb737-de5d-44d5-aaa5-f31246b98b09
# ╠═77862b58-3162-420c-b34e-183d9337b564
# ╠═b44b2a30-d948-44b9-8804-8517eca4ce3c
# ╟─0a6ea318-7414-4943-95ba-b8144984d62d
# ╟─2cf1e8fe-7797-4dc8-85b4-c3a41f1445d6
# ╠═972d3b67-9b8d-46b0-83e1-540b76175516
# ╟─b874c254-222a-4f74-b924-cc785482189b
# ╠═4cd5e33e-964c-46cb-96d5-a868a328bfcf
# ╟─bd78472f-de07-43ad-a21e-85ec6b8a6cce
# ╠═ee2556f0-3c5e-4f85-92ac-37f7422c7057
# ╟─063da746-112c-4c6f-8427-606406870484
# ╠═ee89f54a-823c-4a78-9682-df4cad4349a2
# ╠═a6f1252a-1518-46db-b52d-6c766f913186
# ╟─23c26188-aa54-467b-99a6-be7f40a9655e
# ╠═7a26b34f-bb2d-4834-961b-52aadeb7b0d1
# ╠═d24213c9-9328-421d-b3a1-cc00320c8d80
# ╠═34571dbf-fc63-4a50-8ea5-9ff12a12e1ac
# ╠═ea6d6865-2470-481b-8d1b-4cd0e02468eb
# ╟─306ad617-7e12-4582-a250-f4817004cc5c
# ╠═39ad81dd-391f-46e8-84b6-bfa8456d30b8
# ╟─2a7cee5a-9caa-4760-8923-e9362329f65a
# ╠═80d53fb6-5052-473c-8f17-8ec4a5c88889
