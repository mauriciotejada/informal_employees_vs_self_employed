function initial_param_log_normal(Xmean, Xvar)

    mu = log(Xmean^2 / (sqrt(Xvar + Xmean^2)))
    sig = sqrt(log((1+(Xvar/(Xmean^2)))))

    return mu, sig

end