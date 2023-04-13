
var importObject = {
    env: {
        emscripten_random: function()
        {
            return Math.random();
        },
    },
    wasi_snapshot_preview1: {}
};

WebAssembly.instantiateStreaming(fetch('browser_sim.wasm'), importObject)
.then((results) =>
{
    var getRandomCoordinate = results.instance.exports.getRandomCoordinate;
    var parrotDouble = results.instance.exports.parrotDouble;
    var getTheta = results.instance.exports.getTheta;
    var getThetaDot = results.instance.exports.getThetaDot;
    var getLength = results.instance.exports.getLength;
    var resetAcrobot = results.instance.exports.resetAcrobot;
    var evolveAcrobot = results.instance.exports.evolveAcrobot;
    var getTime = results.instance.exports.getTime;
    var resetAcrobotState = results.instance.exports.resetAcrobotState;
    var applyTorque = results.instance.exports.applyTorque;
    var getTorque = results.instance.exports.getTorque;
    var increaseFriction = results.instance.exports.increaseFriction;
    var getFriction = results.instance.exports.getFriction;
    var getEnergy = results.instance.exports.getEnergy;
    var freezeTorque = results.instance.exports.freezeTorque;
    var holdTorque = results.instance.exports.holdTorque;

    const jointWithFriction = 1;
    var applyJointLock = false;
    var refLockValue = 0.0;
    var applyJointHold = false;
    var refHoldValue = 0.0;
    
    function keyDownEvent(e)
    {
        var code = e.keyCode;
        var key = e.key;

        if (key == 'z' || key == 'Z') {
            resetAcrobotState(0.0, 0.0, 0.0, 0.0);
            applyJointLock = false;
            applyJointHold = false;
        }

        if (key == 'd' || key == 'D') {
            resetAcrobotState(0.01 * (Math.random() - 0.5) - Math.PI / 2, 
                              0.01 * (Math.random() - 0.5) - Math.PI / 2, 
                              0.0, 
                              0.0);
            applyJointLock = false;
            applyJointHold = false;
        }

        if (key == 'u' || key == 'U') {
            resetAcrobotState(0.01 * (Math.random() - 0.5) + Math.PI / 2, 
                              0.01 * (Math.random() - 0.5) + Math.PI / 2, 
                              0.0, 
                              0.0);
            applyJointLock = false;
            applyJointHold = false;
        }

        if (key == 'r' || key == 'R') {
            resetAcrobotState(2.0 * Math.PI * Math.random(), 
                              2.0 * Math.PI * Math.random(), 
                              0.0, 
                              0.0);
            applyJointLock = false;
            applyJointHold = false;
        }

        if ((key == 'f' || key == 'F') && !applyJointHold) {
            applyJointLock = !applyJointLock;
            applyTorque(0.0);
        }

        if ((key == 'h' || key == 'H') && !applyJointLock) {
            applyJointHold = !applyJointHold;
            applyTorque(0.0);
        }
    
        if (code == 39) // right
        {
            increaseFriction(jointWithFriction, 0.01);
        }
        else if (code == 37) // left
        {
            increaseFriction(jointWithFriction, -0.01);
        }
        else if (code == 38) // up
        {
            applyTorque(1.0);
            applyJointLock = false;
            applyJointHold = false;
        }
        else if (code == 40) // down
        {
            applyTorque(-1.0);
            applyJointLock = false;
            applyJointHold = false;
        }
    }

    function keyUpEvent(e)
    {
        var code = e.keyCode;
        if (code == 39) // right
        {

        }
        else if (code == 37) // left
        {

        }
        else if (code == 38) // up
        {
            applyTorque(0.0);
        }
        else if (code == 40) // down
        {
            applyTorque(0.0);
        }
    }

    console.log(parrotDouble(-1.2345));
    console.log(getRandomCoordinate());

    console.log([getLength(0), getLength(1)]);
    
    resetAcrobot();

    console.log([getLength(0), getLength(1)]);

    function draw(canvas, ctx)
    {
        ctx.fillStyle = 'rgb(240, 240, 240)';
        ctx.fillRect(0, 0, canvas.width, canvas.height);

        const L1 = getLength(0);
        const L2 = getLength(1);
        const th1 = getTheta(0);
        const th2 = getTheta(1);

        const cmx1 = L1 * Math.cos(th1) / 2.0;
        const cmy1 = L1 * Math.sin(th1) / 2.0;
        const cmx2 = 2.0 * cmx1 + L2 * Math.cos(th2) / 2.0;
        const cmy2 = 2.0 * cmy1 + L2 * Math.sin(th2) / 2.0;

        const xp0 = canvas.width / 2.0;
        const yp0 = canvas.height / 2.0;
        const ppu = canvas.height / 6.0; // use pixel(x) = xp0 + x * ppu, and pixel(y) = yp0 - y * ppu

        function pixelCoord(x, y) {
          return [xp0 + x * ppu, yp0 - y * ppu];
        }

        function drawArm(cx, cy, L, theta, W) {

          const vx = Math.cos(theta);
          const vy = Math.sin(theta);
          const tx = -vy;
          const ty = vx;

          const hW = W / 2.0;
          const hL = L / 2.0;

          var pix1 = pixelCoord(cx + hL * vx - hW * tx, cy + hL * vy - hW * ty);
          var pix2 = pixelCoord(cx + hL * vx + hW * tx, cy + hL * vy + hW * ty);
          var pix3 = pixelCoord(cx - hL * vx + hW * tx, cy - hL * vy + hW * ty);
          var pix4 = pixelCoord(cx - hL * vx - hW * tx, cy - hL * vy - hW * ty);

          ctx.lineWidth = 3.0;
          ctx.strokeStyle = 'rgb(0, 0, 0)';
          ctx.fillStyle = 'rgb(100, 0, 0)';

          ctx.beginPath();
          ctx.moveTo(pix1[0], pix1[1]);
          ctx.lineTo(pix2[0], pix2[1]);
          ctx.lineTo(pix3[0], pix3[1]);
          ctx.lineTo(pix4[0], pix4[1]);
          ctx.lineTo(pix1[0], pix1[1]);
          ctx.stroke();
          ctx.closePath();
          ctx.fill();
        }

        const W = 0.075;

        drawArm(cmx1, cmy1, L1, th1, W);
        drawArm(cmx2, cmy2, L2, th2, W);

        ctx.fillStyle = 'rgb(0, 0, 0)';
        ctx.font = '16px Arial bold';

        const simTime = getTime();
        ctx.fillText('time: ' + simTime.toFixed(3) + ' [s]', 10.0, 20.0);

        const torque = getTorque();
        ctx.fillText('torque: ' + torque.toFixed(3) + ' [N m]', 10.0, 40.0);

        const mu = getFriction(jointWithFriction);
        ctx.fillText('friction: ' + mu.toFixed(3) + ' [N m s]', 10.0, 60.0);

        const E = getEnergy();
        ctx.fillText('energy: ' + E.toFixed(3) + ' [J]', 10.0, 80.0);

        if (applyJointLock) {
            ctx.fillText('locked at: ' + refLockValue.toFixed(3) + ' [rad]', 10.0, 100.0);
        }

        if (applyJointHold) {
            ctx.fillText('hold at: ' + refHoldValue.toFixed(3) + ' [rad]', 10.0, 100.0);
        }
    }

    const lock_omega_nought = 2.0 * Math.PI / 0.005;
    const lock_feedback_zeta = 1.0;

    const hold_omega_nought = 2.0 * Math.PI / 0.010;
    const hold_feedback_zeta = 1.50;
    
    var startTime = Date.now();
    var delta = 0;
    const dt = 0.001;
    
    function main()
    {
        window.requestAnimationFrame(main);
    
        var currentTime = Date.now();
        var elapsedTime = currentTime - startTime;
        startTime = currentTime;

        var elapsedTimeSeconds = elapsedTime * 1.0e-3;
        delta += elapsedTimeSeconds;

        while (delta >= 0.0) {
            if (applyJointLock) {
                applyTorque(freezeTorque(refLockValue, lock_omega_nought, lock_feedback_zeta));
            } else {
                refLockValue = getTheta(0) - getTheta(1);
            }
            if (applyJointHold) {
                applyTorque(holdTorque(refHoldValue, hold_omega_nought, hold_feedback_zeta));
            } else {
                refHoldValue = getTheta(1);
            }
            evolveAcrobot(dt);
            delta -= dt;
        }
    
        var canvas = document.getElementById('canvas');
        if (canvas.getContext) {
            var ctx = canvas.getContext('2d');
            //ctx.save();
            ctx.clearRect(0, 0, canvas.width, canvas.height);
            draw(canvas, ctx);
            //ctx.restore();
        }
    }
    
    window.addEventListener('keydown', keyDownEvent);
    window.addEventListener('keyup', keyUpEvent);

    window.requestAnimationFrame(main); 

});
