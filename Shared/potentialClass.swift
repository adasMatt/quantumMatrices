//
//  potentialClass.swift
//  quantumMatrices
//
//  Created by Matthew Adas on 3/26/21.
//

import Foundation

class potentialClass: ObservableObject {
    
    var potentialArray: [Double] = []
    var xArray: [Double] = []
    
    func createPotential (potentialType: UInt32, startingX: Double, endingX: Double, stepSize: Double) -> Bool {
        
        //let length = endingX - startingX
        
        switch potentialType {
        case 1:
            for x in stride(from: startingX, through: endingX, by: stepSize) {
                
                xArray.append(x)
                
                let potential = 4.0 * 1.3 * x
                potentialArray.append(potential)
                
            }
        default:
            for x in stride(from: startingX, through: endingX, by: stepSize) {
                
                xArray.append(x)
                
                let potential = 0.0
                potentialArray.append(potential)
                
            }
        }
        
        
        
        if potentialArray.count > 0 {
            return true
        }
        
        else {
            return false
        }
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    func destroyPotential() {
        
        potentialArray = []
        xArray = []
        
    }
    
    
}
