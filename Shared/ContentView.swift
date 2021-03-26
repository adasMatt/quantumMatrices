//
//  ContentView.swift
//  Shared
//
//  Created by Matthew Adas on 3/23/21.
//
// need to fill Eigenvector array instead of returnString
// need E instead of returnString


import SwiftUI
import Accelerate

struct ContentView: View {
    
    @State var resultsString = ""
    @State var numberOfWavefunctionsString = "3"
    @State var startXString = "0.0"
    @State var endXString = "10.0"
    @State var stepSizeString = "0.005"
    @State var potentialTypeString = "1"
    
    @ObservedObject var packUnpack = packUnpackClass()
    @State var wavefunctionArray :[wavefunctions] = []
    @State var currentPotential = potentialClass()
    @State var hamiltonian = hamiltonianClass()
    
    typealias energyAndEigenvectorSolutions = [(E: Double, EigenVector: [Double])]
    
    //let hPlanck = 4.135667696e-15                            // eV*s ...wikipedia
    let electronMass = 0.51099895000                        // ev/c^2 ...wikipedia
    let hbarSquaredOverElectronMass = 7.61996423107385308868  // ev * ang^2
    
    
    var body: some View {
        VStack {

            TextEditor(text: $numberOfWavefunctionsString)
            TextEditor(text: $potentialTypeString)
            TextEditor(text: $startXString)
            TextEditor(text: $endXString)
            TextEditor(text: $stepSizeString)
            TextEditor(text: $resultsString)

            Button("Manipulate Matrices", action: runCalculation)
            
        }
       // .frame(minHeight: 300, maxHeight: 800)
       // .frame(minWidth: 300, maxWidth: 800)
        .padding()
    }
    
    // is this going to be a matrix?
    func createPsi (numberOfWavefunctions: UInt32, startingX: Double, endingX: Double, stepSize: Double) {
        
        for i in 0..<numberOfWavefunctions {
            
            let nthWavefunction = wavefunctions()
            let _ = nthWavefunction.createWaveFunction(nPassed: i+1, startingX: startingX, endingX: endingX, stepSize: stepSize)
            wavefunctionArray.append(nthWavefunction)
            
        }
        
    }
    
    func createV (potentialType: UInt32, startingX: Double, endingX: Double, stepSize: Double) {
        
        currentPotential.destroyPotential()
        let _ = currentPotential.createPotential(potentialType: potentialType, startingX: startingX, endingX: endingX, stepSize: stepSize)
        
    }
    
    func createHam (numberOfWavefunctions: UInt32) {
        hamiltonian.destroyHamiltonian()
        hamiltonian.createHamiltonian(n: UInt32(numberOfWavefunctionsString)!)
    }
    
    func runCalculation () {
        startCalculation(numberOfWavefunctions: UInt32(numberOfWavefunctionsString)!, potentialType: UInt32(potentialTypeString)!, startingX: Double(startXString)!, endingX: Double(endXString)!, stepSize: Double(stepSizeString)!)
    }
    
    func startCalculation (numberOfWavefunctions: UInt32, potentialType: UInt32, startingX: Double, endingX: Double, stepSize: Double) {
        
        /* Calculate Dot Product of 2 Vectors */
        /*
        let q: [Double] = [1, 2, 3]
        let r: [Double] = [0.25, 0.25, 0.25]
        */
        
        let length = endingX - startingX
        
        wavefunctionArray.removeAll()
        currentPotential.destroyPotential()
        
        // calculate currentPotential
        createV(potentialType: potentialType, startingX: startingX, endingX: endingX, stepSize: stepSize)
        var psiByV: [Double] = Array(repeating: 0.0, count: currentPotential.potentialArray.count)
        var result: [Double] = Array(repeating: 0.0, count: currentPotential.potentialArray.count)
        
        // calclate waveFunctionarray?
        createPsi(numberOfWavefunctions: numberOfWavefunctions, startingX: startingX, endingX: endingX, stepSize: stepSize)
        print(wavefunctionArray.count)
        
        // hamiltonian
        createHam(numberOfWavefunctions: numberOfWavefunctions)
        
        for i in 0..<wavefunctionArray.count {
            for j in 0..<wavefunctionArray.count {
                //print(wavefunctionArray[i].psiArray)
                //print(currentPotential.potentialArray)
                
                
                vDSP_vmulD(wavefunctionArray[i].psiArray, 1, currentPotential.potentialArray, 1, &psiByV, 1, vDSP_Length(currentPotential.potentialArray.count))
                
                vDSP_vmulD(wavefunctionArray[j].psiArray, 1, psiByV, 1, &result, 1, vDSP_Length(currentPotential.potentialArray.count))
                
                
                
                var integral = result.mean * length
                
                
                if i == j {
                    integral += wavefunctionArray[i].energy
                    
                }
                
                hamiltonian[i, j] = integral
                
                
                
            }
            
        }
        
        
        // Converts a 2D array into a linear array in FORTRAN Column Major Format (see packingClass file)
        let flatArray: [Double] = packUnpack.pack2dArray(arr: hamiltonian.hamiltonian, rows: Int(wavefunctionArray.count), cols: Int(wavefunctionArray.count))
        
        //print(hamiltonian.hamiltonian)
        // comment these two lines out eventually
        let myString = calculateEigenvaluesString(arrayForDiagonalization: flatArray)
        print(myString)
        
        
        
    }
    /*
    func performMatrixOperations(){
    //func performMatrixOperations(){
        
        //let constants = pow(hPlanck, 2.0) / (8.0 * 1.0) // calculate constants that multiply with momentum squared, unit length 1.0^2.0 as a placeholder for L, 1.0 as unit mass
        
        resultsString += "Add Two Vectors Problem\n\n"
        
        /* Add Two Vectors */
        
        let a: [Double] = [1, 2, 3, 4]
        let b: [Double] = [0.5, 0.25, 0.125, 0.0625]
        var result: [Double] = [0, 0, 0, 0]
        
        resultsString += "["
        resultsString += a.map({"\($0)"}).joined(separator: ", ")
        resultsString += "] + ["
        resultsString += b.map({"\($0)"}).joined(separator: ", ")
        resultsString += "] = "
        
        
        /* vDSP_vaddD(__A: UnsafePointer<Double>, __IA: vDSP_Stride, __B: UnsafePointer<Double>, __IB: vDSP_Stride, __C: UnsafeMutablePointer<Double>, __IC: vDSP_Stride, __N: vDSP_Length) */
        // A: First Vector, IA: Stride of Vector A, B: Second Vector, IB: Stride of Vector B, C: Sum Vector, N: Length of Vectors
        
        vDSP_vaddD(a, 1, b, 1, &result, 1, 4)
        
        resultsString += "["
        resultsString += result.map({"\($0)"}).joined(separator: ", ")
        resultsString += "]\n\n"
        //print(result)
        
        /* Calculate Eigenvalues */
        
        /* Real Eigenvalues */
        
        resultsString += "Real Eigenvalues Problem\n\n"
        
        // need P^2 array and V array?
        //but what about h^2/8mL^2, this is actually just a matrix of n^2
        // L --> xMax? unit length 1.0 for now
        
        // nn member of unperturbed infinite square well energies (electron mass???)
        // nth item in nth array
        
        /*
        /*for array in nSquared {
            item += 1
            
            let newItem = array[item] * constants
            
        }*/
        
        
        // load potential from class? Do I need everything to be 10x10?
        var potentialArray: [Double] = []
        
        //
        var xArray: [Double] = []
        // need a different "realStartingArray" equal to P^2 + V?
        
        //let a = [[[Int]]](repeating:[[Int]](repeating:[Int](repeating:1,count:3), count:3), count:3)
        //var realStartingArray: [[Double]] = [[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]]
        
        // do i need a loop
        /*for array in 0..<9 {
            
            vDSP_vaddD(pSquared2m[array], 1, potentialArray[array], 1, &realStartingArray[array], 1, 10)
            
            vDSP_vadd
        }
        
        //let realStartingArray = pSquared2m + potentialArray       // this made a 20 member array instead of keeping it at 10
        //let realStartingArray = [[2.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0], [4.0, 2.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 4.0, 2.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 4.0, 2.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 4.0, 2.0, 4.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 4.0, 2.0, 4.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 2.0, 4.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 2.0, 4.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 2.0, 4.0], [4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 2.0]]
        
        
     */
        var N = Int32(realStartingArray.count)
        
        // Converts a 2D array into a linear array in FORTRAN Column Major Format (see packingClass file)
        var flatArray: [Double] = packUnpack.pack2dArray(arr: realStartingArray, rows: Int(N), cols: Int(N))
        
        resultsString += "Matrix With Real Eigenvalues\n"
        
        for item in realStartingArray {
        
            resultsString += "\(item)\n"
        }
        
        resultsString += "\n"
        
        // see calculateEigenvalues function for details
        resultsString += calculateEigenvalues(arrayForDiagonalization: flatArray)
        
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
        /* Complex Eigenvalues */
            
        resultsString += "Complex Eigenvalues Problem\n\n"
        
        let complexStartingArray = [[2.0, 4.0], [-4.0, 2.0]]
        
        N = Int32(complexStartingArray.count)
        
        // Converts a 2D array into a linear array in FORTRAN Column Major Format (see packingClass file)
        flatArray = packUnpack.pack2dArray(arr: complexStartingArray, rows: Int(N), cols: Int(N))
        
        resultsString += "Matrix With Complex Eigenvalues\n"
        
        for item in complexStartingArray {
        
            resultsString += "\(item)\n"
        }
        
        resultsString += "\n"
        
        // takes the FORTRAN format of the original array
        resultsString += calculateEigenvalues(arrayForDiagonalization: flatArray)



        
    }
    */
 */
    
    /// calculateEigenvaluesString
    ///
    /// - Parameter arrayForDiagonalization: linear Column Major FORTRAN Array for Diagonalization
    /// - Returns: String consisting of the Eigenvalues and Eigenvectors
    func calculateEigenvaluesString(arrayForDiagonalization: [Double]) -> String {
        /* Integers sent to the FORTRAN routines must be type Int32 instead of Int */
        //var N = Int32(sqrt(Double(startingArray.count)))
        
        var returnString = ""
        
        var N = Int32(sqrt(Double(arrayForDiagonalization.count)))
        var N2 = Int32(sqrt(Double(arrayForDiagonalization.count)))
        var N3 = Int32(sqrt(Double(arrayForDiagonalization.count)))
        var N4 = Int32(sqrt(Double(arrayForDiagonalization.count)))
        
        var flatArray = arrayForDiagonalization
        
        var error : Int32 = 0
        var lwork = Int32(-1)
        // Real parts of eigenvalues
        var wr = [Double](repeating: 0.0, count: Int(N))
        // Imaginary parts of eigenvalues
        var wi = [Double](repeating: 0.0, count: Int(N))
        // Left eigenvectors
        var vl = [Double](repeating: 0.0, count: Int(N*N))
        // Right eigenvectors
        var vr = [Double](repeating: 0.0, count: Int(N*N))
        
        
        /* Eigenvalue Calculation Uses dgeev */
        /*   int dgeev_(char *jobvl, char *jobvr, Int32 *n, Double * a, Int32 *lda, Double *wr, Double *wi, Double *vl,
         Int32 *ldvl, Double *vr, Int32 *ldvr, Double *work, Int32 *lwork, Int32 *info);*/
        
        /* dgeev_(&calculateLeftEigenvectors, &calculateRightEigenvectors, &c1, AT, &c1, WR, WI, VL, &dummySize, VR, &c2, LWork, &lworkSize, &ok)    */
        /* parameters in the order as they appear in the function call: */
        /* order of matrix A, number of right hand sides (b), matrix A, */
        /* leading dimension of A, array records pivoting, */
        /* result vector b on entry, x on exit, leading dimension of b */
        /* return value =0 for success*/
        
        
        
        /* Calculate size of workspace needed for the calculation */
        
        var workspaceQuery: Double = 0.0
        dgeev_(UnsafeMutablePointer(mutating: ("N" as NSString).utf8String), UnsafeMutablePointer(mutating: ("V" as NSString).utf8String), &N, &flatArray, &N2, &wr, &wi, &vl, &N3, &vr, &N4, &workspaceQuery, &lwork, &error)
        
        print("Workspace Query \(workspaceQuery)")
        
        /* size workspace per the results of the query */
        
        var workspace = [Double](repeating: 0.0, count: Int(workspaceQuery))
        lwork = Int32(workspaceQuery)
        
        /* Diagolizing the matrix */
        
        dgeev_(UnsafeMutablePointer(mutating: ("N" as NSString).utf8String), UnsafeMutablePointer(mutating: ("V" as NSString).utf8String), &N, &flatArray, &N2, &wr, &wi, &vl, &N3, &vr, &N4, &workspace, &lwork, &error)
        
        
        if (error == 0)
        {
            for index in 0..<wi.count      /* transform the returned matrices to eigenvalues and eigenvectors */
            {
                if (wi[index]>=0.0)
                {
                    returnString += "Eigenvalue\n\(wr[index]) + \(wi[index])i\n\n"
                }
                else
                {
                    returnString += "Eigenvalue\n\(wr[index]) - \(fabs(wi[index]))i\n\n"
                }
                
                returnString += "Eigenvector\n"
                returnString += "["
                
                
                /* To Save Memory dgeev returns a packed array if complex */
                /* Must Unpack Properly to Get Correct Result
                 
                 VR is DOUBLE PRECISION array, dimension (LDVR,N)
                 If JOBVR = 'V', the right eigenvectors v(j) are stored one
                 after another in the columns of VR, in the same order
                 as their eigenvalues.
                 If JOBVR = 'N', VR is not referenced.
                 If the j-th eigenvalue is real, then v(j) = VR(:,j),
                 the j-th column of VR.
                 If the j-th and (j+1)-st eigenvalues form a complex
                 conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and
                 v(j+1) = VR(:,j) - i*VR(:,j+1). */
                
                for j in 0..<N
                {
                    if(wi[index]==0)
                    {
                        
                        returnString += "\(vr[Int(index)*(Int(N))+Int(j)]) + 0.0i, \n" /* print x */
                        
                    }
                    else if(wi[index]>0)
                    {
                        if(vr[Int(index)*(Int(N))+Int(j)+Int(N)]>=0)
                        {
                            returnString += "\(vr[Int(index)*(Int(N))+Int(j)]) + \(vr[Int(index)*(Int(N))+Int(j)+Int(N)])i, \n"
                        }
                        else
                        {
                            returnString += "\(vr[Int(index)*(Int(N))+Int(j)]) - \(fabs(vr[Int(index)*(Int(N))+Int(j)+Int(N)]))i, \n"
                        }
                    }
                    else
                    {
                        if(vr[Int(index)*(Int(N))+Int(j)]>0)
                        {
                            returnString += "\(vr[Int(index)*(Int(N))+Int(j)-Int(N)]) - \(vr[Int(index)*(Int(N))+Int(j)])i, \n"
                            
                        }
                        else
                        {
                            returnString += "\(vr[Int(index)*(Int(N))+Int(j)-Int(N)]) + \(fabs(vr[Int(index)*(Int(N))+Int(j)]))i, \n"
                            
                        }
                    }
                }
                
                /* Remove the last , in the returned Eigenvector */
                returnString.remove(at: returnString.index(before: returnString.endIndex))
                returnString.remove(at: returnString.index(before: returnString.endIndex))
                returnString.remove(at: returnString.index(before: returnString.endIndex))
                returnString += "]\n\n"
            }
        }
        else {print("An error occurred\n")}
        
        return (returnString)
    }
    
    /// calculateEigenvalues
    ///
    /// - Parameter arrayForDiagonalization: linear Column Major FORTRAN Array for Diagonalization
    /// - Returns: String consisting of the Eigenvalues and Eigenvectors
    
    // typealias energyAndEigenvectorSolutions = [(E: Double, EigenVector: [Double])]
    func calculateEigenvalues(arrayForDiagonalization: [Double]) -> energyAndEigenvectorSolutions {
        /* Integers sent to the FORTRAN routines must be type Int32 instead of Int */
        //var N = Int32(sqrt(Double(startingArray.count)))
        
        var returnString = ""
        
        var N = Int32(sqrt(Double(arrayForDiagonalization.count)))
        var N2 = Int32(sqrt(Double(arrayForDiagonalization.count)))
        var N3 = Int32(sqrt(Double(arrayForDiagonalization.count)))
        var N4 = Int32(sqrt(Double(arrayForDiagonalization.count)))
        
        var flatArray = arrayForDiagonalization
        var solutionsMatrix: energyAndEigenvectorSolutions                      // ??????????????????????????????????????????????????
        
        var error : Int32 = 0
        var lwork = Int32(-1)
        // Real parts of eigenvalues
        var wr = [Double](repeating: 0.0, count: Int(N))
        // Imaginary parts of eigenvalues
        var wi = [Double](repeating: 0.0, count: Int(N))
        // Left eigenvectors
        var vl = [Double](repeating: 0.0, count: Int(N*N))
        // Right eigenvectors
        var vr = [Double](repeating: 0.0, count: Int(N*N))
        
        
        /* Eigenvalue Calculation Uses dgeev */
        /*   int dgeev_(char *jobvl, char *jobvr, Int32 *n, Double * a, Int32 *lda, Double *wr, Double *wi, Double *vl,
         Int32 *ldvl, Double *vr, Int32 *ldvr, Double *work, Int32 *lwork, Int32 *info);*/
        
        /* dgeev_(&calculateLeftEigenvectors, &calculateRightEigenvectors, &c1, AT, &c1, WR, WI, VL, &dummySize, VR, &c2, LWork, &lworkSize, &ok)    */
        /* parameters in the order as they appear in the function call: */
        /* order of matrix A, number of right hand sides (b), matrix A, */
        /* leading dimension of A, array records pivoting, */
        /* result vector b on entry, x on exit, leading dimension of b */
        /* return value =0 for success*/
        
        
        
        /* Calculate size of workspace needed for the calculation */
        
        var workspaceQuery: Double = 0.0
        dgeev_(UnsafeMutablePointer(mutating: ("N" as NSString).utf8String), UnsafeMutablePointer(mutating: ("V" as NSString).utf8String), &N, &flatArray, &N2, &wr, &wi, &vl, &N3, &vr, &N4, &workspaceQuery, &lwork, &error)
        
        print("Workspace Query \(workspaceQuery)")
        
        /* size workspace per the results of the query */
        
        var workspace = [Double](repeating: 0.0, count: Int(workspaceQuery))
        lwork = Int32(workspaceQuery)
        
        /* Diagolizing the matrix */
        
        dgeev_(UnsafeMutablePointer(mutating: ("N" as NSString).utf8String), UnsafeMutablePointer(mutating: ("V" as NSString).utf8String), &N, &flatArray, &N2, &wr, &wi, &vl, &N3, &vr, &N4, &workspace, &lwork, &error)
        
        var eigValueDouble = 0.0
        var eigVectorDoubleArray: [Double] = []
        
        if (error == 0)
        {
            for index in 0..<wi.count      /* transform the returned matrices to eigenvalues and eigenvectors */
            {
                
                // should  I change this entirely to if wi[index] == 0 ... else not zero?
                
                if (wi[index]>=0.0)
                {
                    returnString += "Eigenvalue\n\(wr[index]) + \(wi[index])i\n\n"
                    
                    // need E instead of returnString
                    // typealias energyAndEigenvectorSolutions = [(E: Double, EigenVector: [Double])]
                    eigValueDouble = wr[index]
                }
                else
                {
                    returnString += "Eigenvalue\n\(wr[index]) - \(fabs(wi[index]))i\n\n"
                }
                
                returnString += "Eigenvector\n"
                returnString += "["
                
                
                /* To Save Memory dgeev returns a packed array if complex */
                /* Must Unpack Properly to Get Correct Result
                 
                 VR is DOUBLE PRECISION array, dimension (LDVR,N)
                 If JOBVR = 'V', the right eigenvectors v(j) are stored one
                 after another in the columns of VR, in the same order
                 as their eigenvalues.
                 If JOBVR = 'N', VR is not referenced.
                 If the j-th eigenvalue is real, then v(j) = VR(:,j),
                 the j-th column of VR.
                 If the j-th and (j+1)-st eigenvalues form a complex
                 conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and
                 v(j+1) = VR(:,j) - i*VR(:,j+1). */
                
                for j in 0..<N
                {
                    // need to fill Eigenvector array instead of returnString
                    // typealias energyAndEigenvectorSolutions = [(E: Double, EigenVector: [Double])]
                    if(wi[index]==0)
                    {
                        
                        returnString += "\(vr[Int(index)*(Int(N))+Int(j)]) + 0.0i, \n" /* print x */
                        
                    }
                    else if(wi[index]>0)
                    {
                        if(vr[Int(index)*(Int(N))+Int(j)+Int(N)]>=0)
                        {
                            returnString += "\(vr[Int(index)*(Int(N))+Int(j)]) + \(vr[Int(index)*(Int(N))+Int(j)+Int(N)])i, \n"
                        }
                        else
                        {
                            returnString += "\(vr[Int(index)*(Int(N))+Int(j)]) - \(fabs(vr[Int(index)*(Int(N))+Int(j)+Int(N)]))i, \n"
                        }
                    }
                    else
                    {
                        if(vr[Int(index)*(Int(N))+Int(j)]>0)
                        {
                            returnString += "\(vr[Int(index)*(Int(N))+Int(j)-Int(N)]) - \(vr[Int(index)*(Int(N))+Int(j)])i, \n"
                            
                        }
                        else
                        {
                            returnString += "\(vr[Int(index)*(Int(N))+Int(j)-Int(N)]) + \(fabs(vr[Int(index)*(Int(N))+Int(j)]))i, \n"
                            
                        }
                    }
                }
                
                /* Remove the last , in the returned Eigenvector */
                returnString.remove(at: returnString.index(before: returnString.endIndex))
                returnString.remove(at: returnString.index(before: returnString.endIndex))
                returnString.remove(at: returnString.index(before: returnString.endIndex))
                returnString += "]\n\n"
            }
        }
        else {print("An error occurred\n")}
        
        return (returnString)
    }
    
    
    

}

struct ContentView_Previews: PreviewProvider {
    static var previews: some View {
        ContentView()
    }
}

