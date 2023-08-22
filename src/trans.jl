function trans(x)

    y = similar(x)

    y[1:6] = exp.(x[1:6])
    y[7]  = x[7]
    y[8]  = exp.(x[8])
    y[9]  = x[9]
    y[10] = exp.(x[10])
    y[11] = x[11]
    y[12] = exp.(x[12])
    y[13] = exp.(x[13])
    y[14] = exp.(x[14])

    return y

end

function trans_inv(x)

    y = similar(x)

    y[1:6] = log.(x[1:6])
    y[7]  = x[7]
    y[8]  = log.(x[8])
    y[9]  = x[9]
    y[10] = log.(x[10])
    y[11] = x[11]
    y[12] = log.(x[12])
    y[13] = log.(x[13])
    y[14] = log.(x[14])

    return y

end

function trans_r(x)

    y = similar(x)

    y[1:4] = exp.(x[1:4])
    y[5]  = x[5]
    y[6]  = exp.(x[6])
    y[7]  = x[7]
    y[8] = exp.(x[8])
    y[9] = exp.(x[9])
    y[10] = exp.(x[10])

    return y

end

function trans_inv_r(x)

    y = similar(x)

    y[1:4] = log.(x[1:4])
    y[5]  = x[5]
    y[6]  = log.(x[6])
    y[7]  = x[7]
    y[8] = log.(x[8])
    y[9] = log.(x[9])
    y[10] = log.(x[10])

    return y

end