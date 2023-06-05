#ifndef TMP_TMP_H_
#define TMP_TMP_H_

namespace tmp
{
  int add(int, int);
}

#include <array>
#include <iostream>
#include <memory>
#include <vector>

struct EulerParams
{
  double gamma;
  double ref_temperature;
  double ref_density;
  double ref_pressure;
  class Builder
  {
    double gamma;
    double ref_temperature;
    double ref_density;
    double ref_pressure;

   public:
    Builder& set_gamma(double i_gamma)
    {
      this->gamma = i_gamma;
      return *this;
    }
    Builder& set_ref_temperature(double i_ref_temperature)
    {
      this->ref_temperature = i_ref_temperature;
      return *this;
    }
    Builder& set_ref_density(double i_ref_density)
    {
      this->ref_density = i_ref_density;
      return *this;
    }
    Builder& set_ref_pressure(double i_ref_pressure)
    {
      this->ref_pressure = i_ref_pressure;
      return *this;
    }
    EulerParams build()
    {
      return EulerParams{ gamma, ref_temperature, ref_density, ref_pressure };
    }
  };
};

template<size_t dim, size_t nstate, typename real>
class ConvectivePhysics
{
 public:
  virtual void print() = 0;
  virtual std::array<std::array<real, dim + 2>, dim> compute_flux(const std::array<real, nstate>& state) = 0;
  virtual ~ConvectivePhysics() = default;
};

template<size_t dim, size_t nstate, typename real>
class Euler : public ConvectivePhysics<dim, nstate, real>
{
  const EulerParams params;

 public:
  Euler(const EulerParams& i_params) : params(i_params)
  {
    static_assert(nstate == dim + 2);
  }
  void print() override
  {
    std::cout << "Euler" << std::endl;
    std::cout << "gamma: " << params.gamma << std::endl;
    std::cout << "ref_temperature: " << params.ref_temperature << std::endl;
    std::cout << "ref_density: " << params.ref_density << std::endl;
    std::cout << "ref_pressure: " << params.ref_pressure << std::endl;
  }

  [[nodiscard]] inline std::array<real, dim> compute_velocity(const std::array<real, nstate>& state)
  {
    std::array<real, dim> velocity;
    for (size_t i = 0; i < dim; ++i)
    {
      velocity[i] = state[1 + i] / state[0];
    }
    return velocity;
  }

  [[nodiscard]] inline real compute_energy(const std::array<real, nstate>& state)
  {
    return state[nstate - 1] / (params.gamma - 1) + 0.5 * state[1] * state[1] / state[0];
  }
  [[nodiscard]] inline real compute_pressure(const std::array<real, nstate>& state)
  {
    return (params.gamma - 1) * (state[nstate - 1] - 0.5 * state[1] * state[1] / state[0]);
  }
  [[nodiscard]] inline std::array<std::array<real, nstate>, dim> compute_flux(const std::array<real, nstate>& state) override
  {
    std::array<std::array<real, nstate>, dim> flux;
    const real density = state[0];
    const real energy = compute_energy(state);
    const std::array<real, dim> velocity = compute_velocity(state);
    const real pressure = compute_pressure(state);
    for (size_t idim = 0; idim < dim; ++idim)
    {
      flux[idim][0] = density * velocity[idim];
      for (size_t jdim = 0; jdim < dim; ++jdim)
      {
        flux[idim][1 + jdim] = density * velocity[idim] * velocity[jdim];
      }
      flux[idim][1 + idim] += pressure;
      flux[idim][nstate - 1] = (energy + pressure) * velocity[idim];
    }
    return flux;
  }
};

enum class PhysicsType
{
  EULER,
  NAVIER_STOKES
};

class PhysicsFactory
{
 public:
  template<size_t dim, size_t nstate, typename real>
  static std::unique_ptr<ConvectivePhysics<dim, nstate, real>> create(PhysicsType type)
  {
    switch (type)
    {
      case PhysicsType::EULER:
        return std::make_unique<Euler<dim, nstate, real>>(EulerParams::Builder()
                                                              .set_gamma(1.4)
                                                              .set_ref_temperature(300)
                                                              .set_ref_density(1.225)
                                                              .set_ref_pressure(101325)
                                                              .build());
      default: throw std::runtime_error("Unknown physics type");
    }
  }
};
#endif  // TMP_TMP_H_
