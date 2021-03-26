//
//  hamiltonian.swift
//  quantumMatrices
//
//  Created by Matthew Adas on 3/26/21.
//

import Foundation


class hamiltonianClass: ObservableObject {
    
    var hamiltonian: [[Double]] = [[]]
    
    
    func createHamiltonian(n: UInt32) {
        
        hamiltonian = Array(repeating: Array(repeating: 0, count: Int(n)), count: Int(n))
        
        
        
    }
    
    subscript(row: Int, column: Int) -> Double {
        get {
            // This could validate arguments.
            return hamiltonian[row][column]
        }
        set {
            // This could also validate.
            hamiltonian[row][column] = newValue
        }
    }
    
    // 2D indexable array -- https://stackoverflow.com/questions/25127700/two-dimensional-array-in-swift
    
    //var hamiltonianArray = Array(repeating: Array(repeating: 0, count: 3), count: 3)
    
    func destroyHamiltonian() {
        hamiltonian.removeAll()
    }
    
}
