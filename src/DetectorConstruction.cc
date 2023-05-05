//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software 
// shall cite the following Geant4-DNA collaboration publications:
// Med. Phys. 37 (2010) 4692-4708
// Phys. Med. 31 (2015) 861-874
// The Geant4-DNA web site is available at http://geant4-dna.org
//
/// \file medical/dna/svalue/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4NistManager.hh"
#include "G4Sphere.hh"
#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4RunManager.hh"

#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4UserLimits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
 : G4VUserDetectorConstruction()
{
  //default tracking cut  
  fTrackingCut = 7.4*CLHEP::eV;
  
  // default parameter values
  fWorldRadius = 30*CLHEP::um;
  fCytoThickness = 1*CLHEP::um;
  fNuclRadius = 4*CLHEP::um;
  fECMRadius = 18.09 * CLHEP::um;
  
  DefineMaterials();

  // create commands for interactive definition of the detector  
  fDetectorMessenger = new DetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{ 
  delete fDetectorMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{ 
  G4NistManager* man = G4NistManager::Instance();
  
  fWorldMaterial = man->FindOrBuildMaterial("G4_WATER");
  fCytoMaterial = fNuclMaterial = man->FindOrBuildMaterial("G4_WATER");
  fECMMaterial = man->FindOrBuildMaterial("G4_TISSUE_SOFT_ICRU-4");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
  
G4VPhysicalVolume* DetectorConstruction::Construct()
{
  if(fWorld) return fWorld;
  G4bool checkOverlaps = false;
                
  // tube world
  //
    //World
  G4Box* solidWorld =
      new G4Box("World",                       //its name
          fWorldRadius, fWorldRadius, fWorldRadius);     //its size

  fLogicalWorld =
      new G4LogicalVolume(solidWorld,          //its solid
          fWorldMaterial,           //its material
          "World");            //its name

  fWorld =
      new G4PVPlacement(0,                     //no rotation
          G4ThreeVector(),       //at (0,0,0)
          fLogicalWorld,            //its logical volume
          "World",               //its name
          0,                     //its mother  volume
          false,                 //no boolean operation
          0,                     //copy number
          checkOverlaps);        //overlaps checking                      

  //Spherical ECM
  G4Orb* solidECM =
      new G4Orb("ECM",
          fECMRadius);

  fLogicalECM =
      new G4LogicalVolume(solidECM,            //its solid
          fECMMaterial,             //its material
          "ECM");         //its name

  fECM =
      new G4PVPlacement(0,                       //no rotation
          G4ThreeVector(),         //at (0,0,0)
          fLogicalECM,                //its logical volume
          "ECM",              //its name
          fLogicalWorld,              //its mother  volume
          false,                   //no boolean operation
          0,                       //copy number
          checkOverlaps);          //overlaps checking



  //Here for 13 cells
  // Spherical nucleus
  //
  G4Sphere* 
  sNucl = new G4Sphere("Nucl",                           
                       0., 
                       fNuclRadius, 
                       0., 
                       twopi, 
                       0., 
                       pi);    

  fLogicalNucl = new G4LogicalVolume(sNucl,                        
                                     fNuclMaterial,          
                                    "Nucl");              
                                   
                 

  // Spherical shell for cytoplasm
  //
  G4Sphere* 
  sCyto = new G4Sphere("Cyto",                           
                        fNuclRadius, 
                        fNuclRadius+fCytoThickness, 
                        0., 
                        twopi, 
                        0., 
                        pi);

  fLogicalCyto = new G4LogicalVolume(sCyto,                    
                                     fCytoMaterial,       
                                    "Cyto");       
                                   
 

  //放置细胞
  // 
  G4ThreeVector cell_Vector[13];
  cell_Vector[0] = G4ThreeVector(0 * um, 0 * um, 0 * um);
  cell_Vector[1] = G4ThreeVector(-12.06 * um, 0. * um, 0. * um);
  cell_Vector[2] = G4ThreeVector(12.06 * um, 0. * um, 0. * um);
  cell_Vector[3] = G4ThreeVector(-6.03 * um, 0 * um, -10.44 * um);
  cell_Vector[4] = G4ThreeVector(6.03 * um, 0 * um, -10.44 * um);
  cell_Vector[5] = G4ThreeVector(-6.03 * um, 0 * um, 10.44 * um);
  cell_Vector[6] = G4ThreeVector(6.03 * um, 0 * um, 10.44 * um);
  cell_Vector[7] = G4ThreeVector(0 * um, -10 * um, -5.22 * um);
  cell_Vector[8] = G4ThreeVector(-6.03 * um, -10 * um, 5.22 * um);
  cell_Vector[9] = G4ThreeVector(6.03 * um, -10 * um, 5.22 * um);

  cell_Vector[10] = G4ThreeVector(0 * um, 10 * um, -5.22 * um);
  cell_Vector[11] = G4ThreeVector(-6.03 * um, 10 * um, 5.22 * um);
  cell_Vector[12] = G4ThreeVector(6.03 * um, 10 * um, 5.22 * um);


  /*
按细胞123依次为其所在的位置进行命名，但在logical上不做区分,一定要用copynum
*/
  for (int i = 0; i < 13; i++)
  {

      fNucl = new G4PVPlacement(
          0,
          cell_Vector[i],
          fLogicalNucl,
          "Nucl_",
          fLogicalECM,
          false,
          i,
          checkOverlaps);

      fCyto = new G4PVPlacement(
          0,
          cell_Vector[i],
          fLogicalCyto,
          "Cyto_",
          fLogicalECM,
          false,
          i,
          checkOverlaps);


  }

  //
  
  PrintParameters();
    
  // Tracking cut
  //
  fLogicalNucl->SetUserLimits(new G4UserLimits(DBL_MAX,DBL_MAX,DBL_MAX,
    fTrackingCut));
  fLogicalCyto->SetUserLimits(new G4UserLimits(DBL_MAX,DBL_MAX,DBL_MAX,
    fTrackingCut));
  fLogicalECM->SetUserLimits(new G4UserLimits(DBL_MAX, DBL_MAX, DBL_MAX,
      fTrackingCut));
  fLogicalWorld->SetUserLimits(new G4UserLimits(DBL_MAX,DBL_MAX,DBL_MAX,
    fTrackingCut));
    
  //
  //always return the root volume
  //  
  return fWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::PrintParameters() const
{
  G4cout << "\n---------------------------------------------------------\n";
  G4cout << "---> The tracking cut is set to " 
         << G4BestUnit(fTrackingCut,"Energy") << G4endl;
  G4cout << "---> The World is a sphere of " 
         << G4BestUnit(1000*fNuclRadius,"Length") << "radius of "
         << fWorldMaterial->GetName() << G4endl;
  G4cout << "---> The Nucleus is a sphere of " 
         << G4BestUnit(fNuclRadius,"Length") << "radius of "
         << fWorldMaterial->GetName() << " of mass "
  << G4BestUnit(GetNuclMass(),"Mass") << G4endl;
  G4cout << "---> The Cytoplasm is a spherical shell of thickness " 
         << G4BestUnit(fCytoThickness,"Length") << "of "
         << fWorldMaterial->GetName() << " of mass "
  << G4BestUnit(GetCytoMass(),"Mass") << G4endl;
  G4cout << "\n---------------------------------------------------------\n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetTrackingCut(G4double value)
{
  fTrackingCut = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetNuclRadius(G4double value)
{
  fNuclRadius = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetCytoThickness(G4double value)
{
  fCytoThickness = value;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetWorldMaterial(const G4String& materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial && pttoMaterial != fWorldMaterial) {
    fWorldMaterial = pttoMaterial;
    if(fLogicalWorld) fLogicalWorld->SetMaterial(pttoMaterial);
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetECMMaterial(const G4String& materialChoice)
{
    // search the material by its name   
    G4Material* pttoMaterial =
        G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

    if (pttoMaterial && pttoMaterial != fECMMaterial) {
        fECMMaterial = pttoMaterial;
        if (fLogicalECM) fLogicalECM->SetMaterial(pttoMaterial);
        G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetCytoMaterial(const G4String& materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial && pttoMaterial != fCytoMaterial) {
    fCytoMaterial = pttoMaterial;
    if(fLogicalCyto) fLogicalCyto->SetMaterial(pttoMaterial);
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetNuclMaterial(const G4String& materialChoice)
{
  // search the material by its name   
  G4Material* pttoMaterial =
    G4NistManager::Instance()->FindOrBuildMaterial(materialChoice);

  if (pttoMaterial && pttoMaterial != fNuclMaterial) {
    fNuclMaterial = pttoMaterial;
    if(fLogicalNucl) fLogicalNucl->SetMaterial(pttoMaterial);
    G4RunManager::GetRunManager()->PhysicsHasBeenModified();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
