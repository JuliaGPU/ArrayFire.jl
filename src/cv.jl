### Computer Vision

# Export methods

export  AFFeatures,
        orb,
        sift,
        gloh,
        diffOfGaussians,
        fast,
        harris,
        susan,
        hammingMatcher,
        nearestNeighbour,        
        matchTemplate

# Export constants

export  AF_SAD,
        AF_ZSAD,
        AF_LSAD,
        AF_ZSSD,   
        AF_LSSD,   
        AF_NCC,   
        AF_ZNCC,    
        AF_SHD

# Feature Type

immutable AFFeatures
    ptr::Ptr{Void}
end

# Constants

AF_SAD = 0   
AF_ZSAD = 1
AF_LSAD = 2
AF_ZSSD = 3   
AF_LSSD = 4   
AF_NCC = 5   
AF_ZNCC = 6    
AF_SHD = 7

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

# Feature Detectors

function diffOfGaussians(a::AFArray, radius1::Int, radius2::Int) 
    out = new_ptr()
    af_dog(out, a, radius1, radius)
    AFArray{backend_eltype(out[])}(out[])
end

function fast(a::AFArray; thr = 20., arc_length = 9, non_max = true, 
                feature_ratio = 0.05, edge = 3)
    out = new_ptr()
    af_fast(out, a, thr, Cuint(arc_length), non_max, feature_ratio, Cuint(edge))
    AFFeatures(out[])
end 

function harris(a::AFArray; max_corners = 500, min_response = 1e-5, 
                sigma = 1., block_size = 0, k_thr = 0.04)
    out = new_ptr()
    af_harris(out, a, Cuint(max_corners), min_response, sigma, Cuint(block_size), k_thr)
    AFFeatures(out[])
end 

function susan(a::AFArray; radius = 3, diff_thr = 32.0, 
                geom_thr = 10.0, feature_ratio = 0.05, edge = 3)
    out = new_ptr()
    af_susan(out, a, Cuint(radius), diff_thr, geom_thr, feature_ratio, Cuint(edge))
    AFFeatures(out[])
end 

# Feature Matchers

function hammingMatcher(a::AFArray, train::AFArray; dist_dim = 0, n_dist = 1)
    idx = out_ptr()
    dist = out_ptr()
    af_hamming_matcher(idx, dist, a, train, dist_dim, Cuint(n_dist))
    AFArray{backend_eltype(idx[])}(idx[]) + 1,
    AFArray{backend_eltype(dist[])}(dist[])
end

function nearestNeighbour(a::AFArray, train::AFArray; dist_dim = 9, 
                            n_dist = 1, dist_type = AF_SSD)
    idx = new_ptr()
    dist = new_ptr()
    af_nearest_neighbour(idx, dist, query, train, dist_dim, n_dist, dist_type)
    AFArray{backend_eltype(idx[])}(idx[]) + 1,
    AFArray{backend_eltype(dist[])}(dist[])
end

function matchTemplate(search_img::AFArray, template_img::AFArray; m_typ = AF_SAD)
    out = new_ptr()
    af_match_template(out, search_img, template_img, m_type)
    AFArray{backend_eltype(out[])}(out[])
end
