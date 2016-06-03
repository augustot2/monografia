using PyPlot
P(x,p) = x.^p
x = collect(linspace(-1,1,101));
for p in [0:2:10;]
  temp = 0.0
  for n in 0:p
    temp += P(x,n)
  end
  plot(x,temp,label ="P = $p")
  legend()
  @time sleep(5)
end
