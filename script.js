const fs = require("fs");

function Slay(A, B1, B2, B3, B4, F, C, eps) { // function for chart equation calculation where F - right part of chart equations 
  //  const C = Array(N).fill(0) // area of concentration
    let max = 0

    do {
        for (let i = 1; i < Nx - 1; i++) {
            for (let j = 1; j < Ny - 1; j++) {
                max = 0; //orogin value of deviation
                const m0 = i + j * Nx //calculate the element number,
                const m1 = m0 + 1; //in the center of the pattern chart equation and find
                const m2 = m0 - 1;   //the values of the numbers of elements in the vicinity of the center
                const m3 = m0 + Nx;  //pattern
                const m4 = m0 - Nx;
                const r = C[m0];
                C[m0] = (F[m0] + B1[m0] * C[m1] + B2[m0] * C[m2] + B3[m0] * C[m3] + B4[m4] * C[m0]) / A[m0]
                const module = Math.abs(r - C[m0]); // new value in concentration in this node

                // console.log(`max = ${max}`);
                // console.log(`module = ${module}`);

                if (max < module) {
                    max = module;
                }  // calc of maxvalue of abs of residual between values of val of substance concentration on previous layer and on current    
            }
        }


        tSum += ht;
    } while(max > eps)

   return C;
}



//--------------------------Константы--------------------------
const Nx = 100; //size chart
const Ny = 100; // same
const hx = 0.1; // spatial coordinates steps
const hy = 0.1; //same
const ht = 0.001; //time step
const lt = 1; //time interval
const N = Nx * Ny; // целочисленная переменная, в которой хранится размерность массивов
const eps = 1e-8; //error of calculation of chart equation
const alfa_X = 0; //coefficients that consists of border conditions - *
const alfa_Y = 0; // *
const betta_X = 0; // *
const betta_Y = 0; // *
const sigma = 0.5; // weight of implicit schema

//--------------------------Массивы--------------------------
let C = Array(N).fill(0); //  one-dimensional real array for describing the concentration field
const mu = Array(N).fill(10); //одномерные вещественные массивы для поля коэффициента турбулентного обмена
const f = Array(N).fill(0); //одномерные вещественные массивы для функции, описывающей интенсивность и распределение источников
const u = Array(N).fill(10); //одномерные вещественные массивы для компонентов вектора скорости
const v = Array(N).fill(0); // то же самое
const O = Array(N).fill(0); //‒ одномерные вещественные массивы для функции, описывающей заполненность ячеек
const A = Array(N).fill(0); // одномерные вещественные массивы для коэффициентов сеточных уравнений на текущем временном слое - *
const B1 = Array(N).fill(0);
const B2 = Array(N).fill(0);
const B3 = Array(N).fill(0);
const B4 = Array(N).fill(0);
const B5 = Array(N).fill(0); // одномерные вещественные массивы для коэффициентов сеточных уравнений на предыдущем временном слое
const B6 = Array(N).fill(0);
const B7 = Array(N).fill(0);
const B8 = Array(N).fill(0);
const B9 = Array(N).fill(0);
const F = Array(N).fill(0); //одномерный вещественный массив для правых частей сеточных уравнений

for (let i = 1; i < Nx - 2; i++) { // alg of save 
    for (let j = 1; j < Ny - 2; j++) {

        const m0 = i + j * Nx;
        O[m0] = 1;

    }
}

const i = Nx / 4;
const j = Ny / 4;
const m0 = i + j * Nx;
f[m0] = 1;

let tSum = 0

while (tSum <= lt + ht / 2) {
    for (let i = 1; i < Nx - 1; i++) {
        for (let j = 1; j < Ny - 1; j++) {

            const m0 = i + j * Nx; //nuber to every pair of i j

            const m1 = m0 + 1; // nodes' numbers
            const m2 = m0 - 1; //вычисляют номер элемента, стоящего в центре шаблона сеточного уравнения 
            const m3 = m0 + Nx;//и находят значения номеров элементов, стоящих в окрестности центра  шаблона
            const m4 = m0 - Nx;
            const m24 = m0 - Nx - 1;

            const q1 = (O[m0] + O[m4]) / 2; // fullness of control areas
            const q2 = (O[m2] + O[m24]) / 2;
            const q3 = (O[m0] + O[m2]) / 2;
            const q4 = (O[m4] + O[m24]) / 2;
            const q0 = (q1 + q2) / 2;

            B1[m0] = (-(u[m1] + u[m0]) / (4 * hx) + (mu[m1] + mu[m0]) / (2 * hx * hx)) * q1; //first, the values of the coefficients in the 
            B2[m0] = ((u[m2] + u[m0]) / (4 * hx) + (mu[m2] + mu[m0]) / (2 * hx * hx)) * q2; //neighborhood are calculated
            B3[m0] = ((-v[m3] + v[m0]) / (4 * hy) + (mu[m3] + mu[m0]) / (2 * hy * hy)) * q3; //center of the pattern without taking into 
            B4[m0] = ((v[m4] + v[m0]) / (4 * hy) + (mu[m4] + mu[m0]) / (2 * hy * hy)) * q4; //account the weight of the scheme

            B6[m0] = (1 - sigma) * B1[m0]; //then the coefficients in the locality of the center are calculated
            B7[m0] = (1 - sigma) * B2[m0]; // template on the previous time layer
            B8[m0] = (1 - sigma) * B3[m0];
            B9[m0] = (1 - sigma) * B4[m0];

            B1[m0] = sigma * B1[m0]; // the weight of the scheme is taken into account for the coefficients
            B2[m0] = sigma * B2[m0];
            B3[m0] = sigma * B3[m0];
            B4[m0] = sigma * B4[m0];

            //The coefficients are calculated for the nodes located in the center of the template on:
            A[m0] = q0 / ht + B1[m0] + B2[m0] + B3[m0] + B4[m0] - sigma * (Math.abs(q1 - q2) * (alfa_X / hx) + Math.abs(q3 - q4) * (alfa_Y / hy)) * mu[m0]; //current
            B5[m0] = q0 / ht - B6[m0] - B7[m0] - B8[m0] - B9[m0] + (1 - sigma) * (Math.abs(q1 - q2) * (alfa_X / hx) + Math.abs(q3 - q4) * (alfa_Y / hy)) * mu[m0]; //previous

            F[m0] = q0 * f[m0] - Math.abs(q1 - q2) * mu[m0] * (betta_X / hx) - Math.abs(q3 - q4) * mu[m0] * (betta_Y / hy) + B5[m0] * C[m0] + B6[m0] * C[m1] + B7[m0] * C[m2] + B8[m0] * C[m3] + B9[m0] * C[m4]; // calculation of values of right parts of chart equations
        }
    }

    C = Slay(A, B1, B2, B3, B4, F, C, eps);

    tSum += ht;
}


// for (let i = 0; i < Nx; i++) {
//     for (let j = 0; j < Ny; j++) {
//         fs.appendFileSync("model5.txt", C[i + j * Nx].toString().replace(".", ",") + " ")
//     }
//     fs.appendFileSync("model5.txt", "\n")
// }

const heatmapInstance = h337.create({
    container: document.querySelector(".heatMap")
    });
    
    const points = [];
    
    for (let i = 1; i < Nx - 1; i++) {
    for (let j = 1; j < Ny - 1; j++) {
    const point = {
    x: i * 25,
    y: j * 15,
    value: C[i + j * Nx]
    };
    
    
    points.push(point);
    }
    }
    
    const data = {
    max: Math.max(...C),
    data: points
    };
    
    heatmapInstance.setData(data);




////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// const fs = require("fs");

// function Slay(A, B1, B2, B3, B4, F, C, eps) {
//     let max = 0

//     while(tSum <= lt + ht / 2) {
//         for (let i = 1; i < Nx - 1; i++) {
//             for (let j = 1; j < Ny - 1; j++) {
//                 max = 0
//                 const m0 = i + j * Nx

//                 const m1 = m0 + 1
//                 const m2 = m0 - 1
//                 const m3 = m0 + Nx
//                 const m4 = m0 - Nx

//                 const r = C[m0]

//                 C[m0] = (F[m0] + B1[m0] * C[m1] + B2[m0] * C[m2] + B3[m0] * C[m3] + B4[m4] * C[m0]) / A[m0]
//                 const module = Math.abs(r - C[m0])

//                 if (max < module) max = module
//             }
//         }
        
//         if(max < eps) return C
//     }

//     return C
// }

// const Nx = 100
// const Ny = 100
// const N = Nx * Ny;
// const T = 100;
// const lt = 1000
// const eps = 1e-8
// const alphaX = 0
// const alphaY = 0
// const bettaX = 0
// const bettaY = 0
// const ht = 100;
// const hx = 1;
// const hy = 1
// const sigma = 0.5;

// const mu = Array(N).fill(10)
// const f = Array(N).fill(0)
// const u = Array(N).fill(10)
// const v = Array(N).fill(0)
// const O = Array(N).fill(0)
// const A = Array(N).fill(0)
// const B1 = Array(N).fill(0)
// const B2 = Array(N).fill(0)
// const B3 = Array(N).fill(0)
// const B4 = Array(N).fill(0)
// const B5 = Array(N).fill(0)
// const B6 = Array(N).fill(0)
// const B7 = Array(N).fill(0)
// const B8 = Array(N).fill(0)
// const B9 = Array(N).fill(0)
// let C = Array(N).fill(0)
// const F = Array(N).fill(0)

// for (let i = 1; i < Nx - 2; i++) {
//     for (let j = 1; j < Ny - 2; j++) {
//         const m0 = i + j * Nx
//         O[m0] = 1
//     }
// }

// const i = Nx / 4
// const j = Ny / 4
// f[i + j * Nx] = 1

// let tSum = 0

// while (tSum <= lt + ht / 2) {
//     for (let i = 1; i < Nx - 1; i++) {
//         for (let j = 1; j < Ny - 1; j++) {
//             const m0 = i + j * Nx
//             const m1 = m0 + 1
//             const m2 = m0 - 1
//             const m3 = m0 + Nx
//             const m4 = m0 - Nx
//             const m24 = m0 - Nx - 1

//             const q1 = (O[m0] + O[m4]) / 2
//             const q2 = (O[m2] + O[m24]) / 2
//             const q3 = (O[m0] + O[m2]) / 2
//             const q4 = (O[m4] + O[m24]) / 2
//             const q0 = (q1 + q2) / 2

//             B1[m0] = (-(u[m1] + u[m0]) / (4 * hx) + (mu[m1] + mu[m0]) / (2 * hx * hx)) * q1
//             B2[m0] = ((u[m2] + u[m0]) / (4 * hx) + (mu[m2] + mu[m0]) / (2 * hx * hx)) * q2
//             B3[m0] = ((-v[m3] + v[m0]) / (4 * hy) + (mu[m3] + mu[m0]) / (2 * hy * hy)) * q3
//             B4[m0] = ((v[m4] + v[m0]) / (4 * hy) + (mu[m4] + mu[m0]) / (2 * hy * hy)) * q4

//             B6[m0] = (1 - sigma) * B1[m0]
//             B7[m0] = (1 - sigma) * B2[m0]
//             B8[m0] = (1 - sigma) * B3[m0]
//             B9[m0] = (1 - sigma) * B4[m0]

//             B1[m0] = sigma * B1[m0]
//             B2[m0] = sigma * B2[m0]
//             B3[m0] = sigma * B3[m0]
//             B4[m0] = sigma * B4[m0]

//             A[m0] = q0 / ht + B1[m0] + B2[m0] + B3[m0] + B4[m0] - sigma * (Math.abs(q1 - q2) * (alphaX / hx) + Math.abs(q3 - q4) * (alphaY / hy)) * mu[m0]; //current
//             B5[m0] = q0 / ht - B6[m0] - B7[m0] - B8[m0] - B9[m0] + (1 - sigma) * (Math.abs(q1 - q2) * (alphaX / hx) + Math.abs(q3 - q4) * (alphaY / hy)) * mu[m0]; //previous

//             F[m0] = q0 * f[m0] - Math.abs(q1 - q2) * mu[m0] * (bettaX / hx) - Math.abs(q3 - q4) * mu[m0] * (bettaY / hy) + B5[m0] * C[m0] + B6[m0] * C[m1] + B7[m0] * C[m2] + B8[m0] * C[m3] + B9[m0] * C[m4]; // calculation of values of right parts of chart equations
//         }
//     }

//     C = Slay(A, B1, B2, B3, B4, F, C, eps)

//     tSum += ht
// }

// for (let i = 0; i < Nx; i++) {
//     for (let j = 0; j < Ny; j++) {
//         fs.appendFileSync("model5.txt", C[i + j * Nx].toString().replace(".", ",") + " ")
//     }
//     fs.appendFileSync("model5.txt", "\n")
// }