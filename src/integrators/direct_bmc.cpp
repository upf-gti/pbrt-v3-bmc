
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */


#include "integrators/direct_bmc.h"
#include "interaction.h"
#include "camera.h"
#include "film.h"
#include "paramset.h"
#include "../GPRender/src/cov_kernels.h"

// GTI

namespace pbrt {

Vector3f random(Float min, Float max) {
    Float x = min + (max - min) * (rand() / (RAND_MAX + 1.0));
    Float y = min + (max - min) * (rand() / (RAND_MAX + 1.0));
    Float z = min + (max - min) * (rand() / (RAND_MAX + 1.0));
    return Vector3f(x, y, z);
};

// Get a random unit vector
Vector3f random_unit_vector() {
    while (true) {
        Vector3f p = random(-1.0, 1.0);
        Float len_sq = p.LengthSquared();
        if (1e-7 < len_sq && len_sq <= 1.0) {
            return p / std::sqrt(len_sq);
        }
    }
}

// Get random vector on an hemishpere facing the Z axis
Vector3f random_on_hemisphere(const void *data) {
    const Vector3f *normal = static_cast<const Vector3f *>(data);
    Vector3f on_unit_sphere = random_unit_vector();
    if (Dot(on_unit_sphere, *normal) >
        0.0) {  // In the same hemisphere as the normal
        return on_unit_sphere;
    } else {
        return -on_unit_sphere;
    }
}

struct sSamplingParams {
    Vector3f normal;
};

// DirectIntegrator Method Definitions
Spectrum DirectBMCIntegrator::Li(const RayDifferential &ray, const Scene &scene,
                               Sampler &sampler, MemoryArena &arena,
                               int depth) const {
    
    SurfaceInteraction isect, random_si;
    bool foundIntersection = scene.Intersect(ray, &isect);
    Spectrum L(0.0), f(0.0);

    // Initialize common variables for Direct integrator
    Vector3f wo = isect.wo;

    isect.ComputeScatteringFunctions(ray, arena);
    if (!isect.bsdf) return isect.Le(wo);

    // Add contribution of each light source
    Vector3f wiLocal, wiWorld;
    Float pdf = 1.0 / (2.0 * PI);
    VisibilityTester visibility;

    // Random angle to rotate the GP directions
    Float alpha = 2.0 * PI * sampler.Get1D();
    Vector3f normal = {0.0, 0.0, 1.0};
    sSamplingParams *sampling_params = new sSamplingParams();
    sampling_params->normal = normal;

    for (uint32_t i = 0; i < num_shading_samples; ++i) {
        wiLocal = random_on_hemisphere(sampling_params);
        wiWorld = Normalize(isect.bsdf->LocalToWorld(wiLocal));

        Ray nextRay = isect.SpawnRay(wiWorld);
        bool foundRandomIntersection = scene.Intersect(nextRay, &random_si);
        if (!foundRandomIntersection) continue;

        f = isect.bsdf->f(wo, wiWorld) * Dot(wiWorld, isect.shading.n);
        L += random_si.Le(-nextRay.d) * f;
    }
        
    L /= (num_shading_samples * pdf);

    
    // Compute emitted light if ray hit an area light source
    L += isect.Le(wo);

    return L;
}

std::unique_ptr<DirectBMCIntegrator> DirectBMCIntegrator::CreateDirectBMCIntegrator(
    const ParamSet &params, std::shared_ptr<Sampler> sampler,
    std::shared_ptr<const Camera> camera) {
    int numSamples = params.FindOneInt("numsamples", 32);
    int np;
    const int *pb = params.FindInt("pixelbounds", &np);
    Bounds2i pixelBounds = camera->film->GetSampleBounds();
    if (pb) {
        if (np != 4)
            Error("Expected four values for \"pixelbounds\" parameter. Got %d.",
                  np);
        else {
            pixelBounds = Intersect(pixelBounds,
                                    Bounds2i{{pb[0], pb[2]}, {pb[1], pb[3]}});
            if (pixelBounds.Area() == 0)
                Error("Degenerate \"pixelbounds\" specified.");
        }
    }

    std::unique_ptr<DirectBMCIntegrator> bmc_integrator =
        std::make_unique<DirectBMCIntegrator>(numSamples, camera, sampler,
                                              pixelBounds);

    bmc_integrator->bmc_list.resize(bmc_integrator->num_bmcs);

    for (uint32_t i = 0; i < bmc_integrator->num_bmcs; ++i) {
        pbrt_kernel::sSobolevParams *sobolev_params =
            new pbrt_kernel::sSobolevParams();
        sobolev_params->s = 1.5f;

        GaussianProcess<Vector3f, SampledSpectrum>::sKernelInfo kernel_info;

        kernel_info.kernel = pbrt_kernel::sobolev;
        kernel_info.kernel_params = sobolev_params;

        GaussianProcess<Vector3f, SampledSpectrum> *gaussian_process =
            new GaussianProcess<Vector3f, SampledSpectrum>(kernel_info, 0.01);

        // Set x number of samples (observation/training points), in our case
        // directions
        std::vector<Vector3f> sample_directions;
        sample_directions.reserve(bmc_integrator->num_shading_samples);

        // Victor's birth year plus offset :)
        srand(1998 + i);

        // Generate x random directions in sphere and store in array
        Vector3f normal = Vector3f(0.0, 0.0, 1.0);
        for (uint32_t s_idx = 0; s_idx < bmc_integrator->num_shading_samples;
             s_idx++) {
            sample_directions.push_back(random_on_hemisphere(&normal));
        }

        // Fill the GP instance with the array of directions (observation
        // points)
        gaussian_process->set_observations(sample_directions, {});

        sSamplingParams *sampling_params = new sSamplingParams();
        sampling_params->normal = normal;
        BMC<Vector3f, SampledSpectrum>::sSamplingInfo sampling_info;
        sampling_info.sample = random_on_hemisphere;
        sampling_info.sample_params = sampling_params;

        bmc_integrator->bmc_list[i] =
            new BMC<Vector3f, SampledSpectrum>(sampling_info, gaussian_process);
    }

    printf("\nExecuting Direct BMC Integrator\n");
    return bmc_integrator;
}

}  // namespace pbrt