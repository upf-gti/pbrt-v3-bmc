
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


#include "integrators/area.h"
#include "interaction.h"
#include "camera.h"
#include "film.h"
#include "paramset.h"
#include "lights/diffuse.h"

// GTI

namespace pbrt {

struct sSamplingParams_area {
    Vector3f corner;
    Vector3f diagonal;
};

inline Vector3f random_on_area(const void *data) {
    const sSamplingParams_area *params = static_cast<const sSamplingParams_area *>(data);
    Vector3f corner = params->corner;
    Vector3f diagonal = params->diagonal;

    Float r1 = rand() / (RAND_MAX + 1.0);
    Float r2 = rand() / (RAND_MAX + 1.0);

    return Vector3f(corner.x + r1 * diagonal.x, corner.y,
                    corner.z + r2 * diagonal.z);
}


// AreaIntegrator Method Definitions
Spectrum AreaIntegrator::Li(const RayDifferential &ray, const Scene &scene,
                               Sampler &sampler, MemoryArena &arena,
                               int depth) const {

    SurfaceInteraction isect, random_si;
    if (!scene.Intersect(ray, &isect)) {
        return Spectrum(0.0);
    }
    
    // Initialize common variables for Direct integrator
    Vector3f wo = isect.wo;
    Spectrum Le = isect.Le(wo);
    isect.ComputeScatteringFunctions(ray, arena);
    //if (!isect.bsdf) return Le;
    if (Le.MaxComponentValue() > 0.0) return Le;
                

    // Initialize common variables for Area integrator
    Spectrum L(0.0), f(0.0);
    Normal3f normal = isect.shading.n;
    sSamplingParams_area *samplingParams = new sSamplingParams_area();
    std::shared_ptr<Light> light = scene.lights[0];
    
    Vector3f wi;
    Float pdf;
    VisibilityTester visibility;
    Le = light->Sample_Li(isect, sampler.Get2D(), &wi, &pdf, &visibility);
    Normal3f light_normal = visibility.P1().n;
    auto *area_light =
        dynamic_cast<const DiffuseAreaLight *>(light.get());
    Bounds3f b = area_light->shape->WorldBound();
    samplingParams->corner = Vector3f(b.pMin.x, b.pMin.y, b.pMin.z);
    samplingParams->diagonal =
        Vector3f(b.pMax.x - b.pMin.x, b.pMax.y - b.pMin.y, b.pMax.z - b.pMin.z);
    Float area = samplingParams->diagonal.x * samplingParams->diagonal.y +
                 samplingParams->diagonal.x * samplingParams->diagonal.z +
                 samplingParams->diagonal.y * samplingParams->diagonal.z;
    Float pdf_area = 1.0 / area;
    Point3f x = isect.p;
    for (int i = 0; i < num_shading_samples; ++i) {

        Vector3f random_point = random_on_area(samplingParams);
        Vector3f vector_to_light = random_point - Vector3f(x.x, x.y, x.z);
        Vector3f wi = Normalize(vector_to_light);
        Float maxT = vector_to_light.Length();
        Ray nextRay = isect.SpawnRay(wi);

        if (!scene.Intersect(nextRay, &random_si)) {
            if (random_si.Le(wi).MaxComponentValue() > 0.0) {
                f = isect.bsdf->f(wo, wi) * Dot(wi, normal);
                Float geo = Dot(-wi, normal) / vector_to_light.LengthSquared();

                L += Le * f * geo / pdf_area;
            }

        }

    }
        
    L /= num_shading_samples;
        
    return L;
}

AreaIntegrator *CreateAreaIntegrator(
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
    return new AreaIntegrator(numSamples, camera, sampler, pixelBounds);
}

}  // namespace pbrt
