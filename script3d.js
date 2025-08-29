const canvas = document.getElementById("my-canvas");
const ctx = canvas.getContext("2d");

const n = 1000;
const dt = 0.02;
const frictionHalfLife = 0.040;
const rMax = 0.4;
const m = 6;
const matrix = makeRandomMatrix();
const forceFactor = 10;

const frictionFactor = Math.pow(0.5, dt / frictionHalfLife);

function makeRandomMatrix() {
    const rows = [];
    for (let i = 0; i < m; i++) {
        const row = [];
        for (let j = 0; j < m; j++) {
            row.push(Math.random() * 2 - 1);
        }
        rows.push(row);
    }
    return rows;
}

const colors = new Int32Array(n);
const positionsX = new Float32Array(n);
const positionsY = new Float32Array(n);
const positionsZ = new Float32Array(n);
const velocitiesX = new Float32Array(n);
const velocitiesY = new Float32Array(n);
const velocitiesZ = new Float32Array(n);

for(let i = 0; i < n; i++) {
    colors[i] = Math.floor(Math.random() * m);
    positionsX[i] = Math.random() * 2 - 1;
    positionsY[i] = Math.random() * 2 - 1;
    positionsZ[i] = Math.random() * 2 - 1;
    velocitiesX[i] = 0;
    velocitiesY[i] = 0;
    velocitiesZ[i] = 0;
}

function updateParticles() {
    // update velocities
    for (let i = 0; i < n; i++) {
        let totalForceX = 0;
        let totalForceY = 0;
        let totalForceZ = 0;
        
        for (let j = 0; j < n; j++) {
            if (j == i) continue;
            const rx = positionsX[j] - positionsX[i];
            const ry = positionsY[j] - positionsY[i];
            const rz = positionsZ[j] - positionsZ[i];
            const r = Math.sqrt(rx*rx + ry*ry + rz*rz);
            if (r > 0 && r < rMax) {
                const f = force(r / rMax, matrix[colors[i]][colors[j]]);
                totalForceX += rx / r * f;
                totalForceY += ry / r * f;
                totalForceZ += rz / r * f;
            }
        }
        
        totalForceX *= rMax * forceFactor;
        totalForceY *= rMax * forceFactor;
        totalForceZ *= rMax * forceFactor;

        velocitiesX[i] *= frictionFactor;
        velocitiesY[i] *= frictionFactor;
        velocitiesZ[i] *= frictionFactor;

        velocitiesX[i] += totalForceX * dt;
        velocitiesY[i] += totalForceY * dt;
        velocitiesZ[i] += totalForceZ * dt;
    }

    // update positions    
    for (let i = 0; i < n; i++) {
        positionsX[i] += velocitiesX[i] * dt;
        positionsY[i] += velocitiesY[i] * dt;
        positionsZ[i] += velocitiesZ[i] * dt;
    }
}

function loop() {
    updateParticles();
    
    ctx.fillStyle = "black";
    ctx.fillRect(0, 0, canvas.width, canvas.height);    
    
    for (let i = 0; i < n; i++) {
        ctx.beginPath();
        const f = 1 / (positionsZ[i] + 2);
        const screenX = (f * positionsX[i] + 1) * 0.5 * canvas.width;
        const screenY = (f * positionsY[i] + 1) * 0.5 * canvas.height;
        ctx.fillStyle = `hsl(${360 * (colors[i] / m)}, 100%, 50%)`; // Corretto qui
        ctx.arc(screenX, screenY, 1, 0, 2 * Math.PI);
        ctx.fill();
    }
    
    requestAnimationFrame(loop);
}

function force(r, a) {
    const beta = 0.3;
    if (r < beta) {
        return r / beta - 1;
    } else if (beta < r && r < 1) {
        return a * (1 - Math.abs(2 * r - 1 - beta) / (1 - beta));
    } else {
        return 0;
    }
}

requestAnimationFrame(loop);
