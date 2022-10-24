/**
 * @file PrototypeModel.h
 * @brief Class providing a prototype model for dark brem
 * @author Michael Revering, University of Minnesota
 * @author Tom Eichlersmith, University of Minnesota
 */

#ifndef G4DARKBREM_PROTOTYPEMODEL_H
#define G4DARKBREM_PROTOTYPEMODEL_H

#include "G4DarkBreM/Parameters.h"

namespace g4db {

/**
 * @class G4DarkBremmsstrahlungModel
 * Abstract class representing a model for dark brem.
 *
 * The model is what actually determines two important things:
 *  1. How the cross section is calculated
 *  2. What the particle change is when the process happens
 *
 * This class is the base class that shows what is necessary
 * for the model to function properly.
 */
class PrototypeModel {
 public:
  /**
   * Constructor
   *
   * Configures the model based on the passed parameters
   *
   * Names the logger after the name for this model.
   */
  PrototypeModel(const framework::config::Parameters& p, bool muons) {
    muons_ = muons;
  }

  /// Destructor, nothing on purpose
  virtual ~PrototypeModel() = default;

  /**
   * Print the configuration of this model
   *
   * Helpful for debugging and keeping the process compliant
   * with the other Geant4 processes.
   */
  virtual void PrintInfo() const = 0;

  /**
   * Calculate the cross section given the input parameters
   *
   * @see G4DarkBremmstrahlung::GetMeanFreePath
   * @param[in] electronKE current electron kinetic energy
   * @param[in] atomicA atomic-mass number for the element the electron is in
   * @param[in] atomicZ atomic-number for the element the electron is in
   * @returns cross section with units incorporated as a G4double
   */
  virtual G4double ComputeCrossSectionPerAtom(G4double electronKE,
                                              G4double atomicA,
                                              G4double atomicZ) = 0;

  /**
   * Generate the change in the particle now that we can assume the interaction
   * is occuring
   *
   * @note The input particleChange has already been cleared and then
   * initialized, so there is no need for the model to do those steps.
   *
   * @see G4DarkBremmstrahlung::PostStepDoIt
   * @param[in,out] particleChange particle change class that stores information
   * @param[in] track current track that needs the change
   * @param[in] step current step of the track
   */
  virtual void GenerateChange(G4ParticleChange& particleChange,
                              const G4Track& track, const G4Step& step) = 0;

 protected:
  /// whether muons (true) or electrons (false) are dark bremming
  bool muons_;

 public:
  /**
   * Factory to create models from a registered name
   */
  class Factory {
   public:
    /// the handle of our models
    using ModelPtr = std::shared_ptr<PrototypeModel>;
    /// signature of a function to dynamically create models
    using ModelMaker = ModelPtr(*)(const Parameters& p, bool muons);
   public:
    /**
     * get the factory instance
     *
     * Using a static function variable gaurantees that the factory
     * is created as soon as it is needed and that it is deleted
     * before the program completes.
     *
     * @returns reference to single Factory instance
     */
    static Factory& get() {
      static Factory the_factory;
      return the_factory;
    }
  
    /**
     * register a new object to be constructible
     *
     * We insert the new object into the library after
     * checking that it hasn't been defined before.
     *
     * @tparam DerivedType object type to declare
     * @return value to define a static variable to force running this function
     *  at library load time. It relates to variables so that it cannot be
     *  optimized away.
     */
    template<typename DerivedType>
    uint64_t declare(const std::string& name) {
      if (library_.find(name) != library_.end()) {
        throw std::runtime_error("Double Declaration: Another model has been named '"+name+"'.");
      }
      library_[name] = &maker<DerivedType>;
      return reinterpret_cast<std::uintptr_t>(&library_);
    }

    /**
     * make a new object by name
     *
     * We look through the library to find the requested object.
     * If found, we create one and return a pointer to the newly
     * created object. If not found, we raise an exception.
     *
     * @throws Exception if the input object name could not be found
     *
     * The arguments to the maker are determined at compiletime
     * using the template parameters of Factory.
     *
     * @param[in] full_name name of class to create, same name as passed to declare
     * @param[in] maker_args parameter pack of arguments to pass on to maker
     *
     * @returns a pointer to the parent class that the objects derive from.
     */
    ModelPtr make(const std::string& full_name,
                      PrototypeConstructorArgs... maker_args) {
      auto lib_it{library_.find(full_name)};
      if (lib_it == library_.end()) {
        throw Exception("Factory","An object named " + full_name +
                         " has not been declared.",false);
      }
      return lib_it->second(maker_args...);
    }
  
    /// delete the copy constructor
    Factory(Factory const&) = delete;
  
    /// delete the assignment operator
    void operator=(Factory const&) = delete;
  
   private:
    /**
     * make a new DerivedType returning a ModelPtr
     *
     * Basically a copy of what 
     * [`std::make_unique`](https://en.cppreference.com/w/cpp/memory/unique_ptr/make_unique) 
     * or 
     * [`std::make_shared`](https://en.cppreference.com/w/cpp/memory/shared_ptr/make_shared)
     * do but with the following changes:
     *  1. constructor arguments defined by the Factory and not here
     *  2. return type is a base pointer and not a derived pointer
     *
     * This is where we required that ModelPtr has the same
     * behavior as STL smart pointers. The ModelPtr class must
     * be able to be constructed from a pointer to a derived class
     * and must take ownership of the new object.
     *
     * @tparam DerivedType type of derived object we should create
     * @param[in] args constructor arguments for derived type construction
     */
    template <typename DerivedType>
    static ModelPtr maker(const Parameters& p, bool muon) {
      return ModelPtr(new DerivedType(p, muon));
    }
  
    /// private constructor to prevent creation
    Factory() = default;
  
    /// library of possible objects to create
    std::unordered_map<std::string, ModelMaker> library_;
  };  // Factory
};  // PrototypeModel

} // namespace g4db

#endif
