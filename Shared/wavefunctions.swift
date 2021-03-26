//
//  wavefunctions.swift
//  quantumMatrices
//
//  Created by Matthew Adas on 3/26/21.
//

import Foundation

class wavefunctions: ObservableObject {
    
    let hbarSquaredOverElectronMass = 7.61996423107385308868  // ev * ang^2
    
    var psiArray: [Double] = []
    var xArray: [Double] = []
    
    var energy: Double = 0.0
    var n: UInt32 = 0
    
    func createWaveFunction (nPassed: UInt32, startingX: Double, endingX: Double, stepSize: Double) -> Bool {
        
        let length = endingX - startingX
        n = nPassed
        energy = hbarSquaredOverElectronMass * pow(Double.pi, 2.0) * pow(Double(n), 2.0) / (2.0 * length * length)
        
        for x in stride(from: startingX, through: endingX, by: stepSize) {
            
            xArray.append(x)
            
            let psi = sqrt(2.0 / length) * sin(Double(n) * Double.pi * x / length)
            psiArray.append(psi)
            
        }
        
        
        if psiArray.count > 0 {
            return true
        }
        
        else {
            return false
        }
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    func destroyWaveFunction() {
        
        psiArray = []
        xArray = []
        energy = 0.0
        
    }
    
    
}
