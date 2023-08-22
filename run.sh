#!/bin/bash

julia --project=@. main.jl

julia --project=@. main_r.jl

julia --threads 8 --project=@. main_boot.jl

julia --threads 8 --project=@. main_boot_r.jl

julia --threads 8 --project=@. main_getb_boot.jl

julia --threads 8 --project=@. main_getb_boot_r.jl
