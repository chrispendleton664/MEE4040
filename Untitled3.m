% newmark method
dadtnxt = dadt + ((1-alpha)*ddadtt+alpha * ddadttnxt) * dt
anxt = a + dadt*dt + ((0.5-beta)*ddadtt + beta*ddadttnxt)*dt^2
