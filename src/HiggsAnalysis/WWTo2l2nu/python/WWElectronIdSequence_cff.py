import FWCore.ParameterSet.Config as cms

from RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi import *

import RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi
egammaIDCutsLoose = RecoEgamma.ElectronIdentification.electronIdCutBasedExt_cfi.eidCutBasedExt.clone()
egammaIDCutsLoose.electronQuality = 'loose'

#from RecoEgamma.ElectronIdentification.electronIdCutBasedClassesExt_cfi import *

#import RecoEgamma.ElectronIdentification.electronIdCutBasedClassesExt_cfi
#egammaIDCutsLoose = RecoEgamma.ElectronIdentification.electronIdCutBasedClassesExt_cfi.eidCutBasedClassesExt.clone()

WWElectronIdSequence = cms.Sequence( egammaIDCutsLoose )


#from RecoEgamma.ElectronIdentification.electronIdCutBasedClassesExt_cfi import *

#import RecoEgamma.ElectronIdentification.electronIdCutBasedClassesExt_cfi
#egammaIDCutsLoose = RecoEgamma.ElectronIdentification.electronIdCutBasedClassesExt_cfi.eidCutBasedClassesExt.clone()

#from RecoEgamma.ElectronIdentification.electronIdLikelihoodExt_cfi import *

#import RecoEgamma.ElectronIdentification.electronIdLikelihoodExt_cfi
#egammaIDLikelihood = RecoEgamma.ElectronIdentification.electronIdLikelihoodExt_cfi.eidLikelihoodExt.clone()

#WWElectronIdSequence = cms.Sequence( egammaIDCutsLoose + egammaIDLikelihood )

#################################
