# all conductors giving cyclotomic fields with degree <= B

using Hecke

if isassigned(ARGS, 1) B = parse(Int, ARGS[1]) end

function conductors_up_to(B)
  b = Int(floor(B/2))
  c = []
  for i in 1:b
    e = euler_phi_inv(2*i)
    if length(e) != 0
      append!(c, e)
    end
  end
  println(length(c))
  println(maximum(c))
  println(c)
  println(findfirst(x->x==120, c))
end

conductors_up_to(B)
