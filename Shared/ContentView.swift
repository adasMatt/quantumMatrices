//
//  ContentView.swift
//  Shared
//
//  Created by Matthew Adas on 3/23/21.
//

//search "multiply psi*(H+V_embed)*psi instead"

import SwiftUI
import Accelerate

struct ContentView: View {
    
    @State var resultsString = ""
    @ObservedObject var packUnpack = packUnpackClass()
    
    var body: some View {
        VStack {


            TextEditor(text: $resultsString)

            Button("Manipulate Matrices", action: performMatrixOperations)
                     
             
            
        }
       // .frame(minHeight: 300, maxHeight: 800)
       // .frame(minWidth: 300, maxWidth: 800)
        .padding()
    }
    
    func performMatrixOperations(){
        
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
        
        
        //do I not need this at all but rather need to make a multiplication function????
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
        
        //change this to multiply psi*(H+V_embed)*psi instead of having numbers?
        //let a = [[[Int]]](repeating:[[Int]](repeating:[Int](repeating:1,count:3), count:3), count:3)
        let realStartingArray = [[2.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0], [4.0, 2.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 4.0, 2.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 4.0, 2.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 4.0, 2.0, 4.0, 0.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 4.0, 2.0, 4.0, 0.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 2.0, 4.0, 0.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 2.0, 4.0, 0.0], [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 2.0, 4.0], [4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 2.0]]
        
        var N = Int32(realStartingArray.count)
        
        var flatArray: [Double] = packUnpack.pack2dArray(arr: realStartingArray, rows: Int(N), cols: Int(N))
        
        resultsString += "Matrix With Real Eigenvalues\n"
        
        for item in realStartingArray {
        
            resultsString += "\(item)\n"
        }
        
        resultsString += "\n"
        
        resultsString += calculateEigenvalues(arrayForDiagonalization: flatArray)
        
        
        /* Complex Eigenvalues */
            
        resultsString += "Complex Eigenvalues Problem\n\n"
        
        let complexStartingArray = [[2.0, 4.0], [-4.0, 2.0]]
        
        N = Int32(complexStartingArray.count)
        
        flatArray = packUnpack.pack2dArray(arr: complexStartingArray, rows: Int(N), cols: Int(N))
        
        resultsString += "Matrix With Complex Eigenvalues\n"
        
        for item in complexStartingArray {
        
            resultsString += "\(item)\n"
        }
        
        resultsString += "\n"
        
        resultsString += calculateEigenvalues(arrayForDiagonalization: flatArray)


        resultsString += "Calculate Dot Product of Two Vectors Problem\n\n"
        
        
        
        /* Calculate Dot Product of 2 Vectors */
        
        let q: [Double] = [1, 2, 3]
        let r: [Double] = [0.25, 0.25, 0.25]
        
        var dotProductScalar = 0.0
        
        /* vDSP_dotprD(__A: UnsafePointer<Double>, __IA: vDSP_Stride, __B: UnsafePointer<Double>, __IB: vDSP_Stride, __C: UnsafeMutablePointer<Double>, __N: vDSP_Length) */
        // A: First Vector, IA: Stride of Vector A, B: Second Vector, IB: Stride of Vector B, C: dotProduct N: Length of Vectors
        
        vDSP_dotprD(q, 1, r, 1, &dotProductScalar, vDSP_Length(q.count))
        
        resultsString += "["
        resultsString += q.map({"\($0)"}).joined(separator: ", ")
        resultsString += "] • ["
        resultsString += r.map({"\($0)"}).joined(separator: ", ")
        resultsString += "] = "
        
        resultsString += "\(dotProductScalar)"
        resultsString += "\n\n"

        
    }
    
    /// calculateEigenvalues
    ///
    /// - Parameter arrayForDiagonalization: linear Column Major FORTRAN Array for Diagonalization
    /// - Returns: String consisting of the Eigenvalues and Eigenvectors
    func calculateEigenvalues(arrayForDiagonalization: [Double]) -> String {
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
        
        /* Calculate the size of the workspace */
        
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
    
    
    
}


struct ContentView_Previews: PreviewProvider {
    static var previews: some View {
        ContentView()
    }
}
