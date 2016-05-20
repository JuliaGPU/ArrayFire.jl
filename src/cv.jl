### Computer Vision

# Export methods

export  AFFeatures,
        orb,
        sift,
        gloh

# Feature Type

immutable AFFeatures
    ptr::Ptr{Void}
end

# Feature Descriptors

function orb(a::AFArray; fast_thr = 20., max_feat = 400, scl_fctr = 1.5, levels = 4, blur_img = false)
    feat = new_ptr()
    desc = new_ptr()
    af_orb(feat, desc, a, fast_thr, Cuint(max_feat), scl_fctr, Cuint(levels), blur_img)
    AFFeatures(feat[]),
    AFArray{backend_eltype(desc[])}(desc[])
end

for (op,fn) in ((:sift, :af_sift), (:gloh, :af_gloh))

    @eval function ($op)(a::AFArray; n_layers = 3, constant_thr = 0.04, edge_thr = 0.04, 
                    init_sigma = 1.6, double_input = true, intensity_scale = 0.00390625, 
                    feature_ratio = 0.05)
        feat = new_ptr()
        desc = new_ptr()
        eval($fn)(feat, desc, a, Cuint(n_layers), contrast_thr, edge_thr, init_sigma, 
                double_input, intensity_scale, feature_ratio)
        AFFeatures(feat[]),
        AFArray{backend_eltype(desc[])}(desc[])
    end

end
