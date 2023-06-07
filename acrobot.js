
var importObject = { env: { }, };

WebAssembly.instantiateStreaming(fetch('browser_sim.wasm'), importObject)
.then((results) =>
{
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
    var setFriction = results.instance.exports.setFriction;
    var getEnergy = results.instance.exports.getEnergy;
    var freezeTorque = results.instance.exports.freezeTorque;
    var holdTorque = results.instance.exports.holdTorque;
    var cosineAlpha = results.instance.exports.cosineAlpha;
    var directionIndicator = results.instance.exports.directionIndicator;

    const refAngleDelta = 2.0 * Math.PI / 24.0;
    const jointWithFriction = 1;
    const jointWithBrake = 0;

    var applyJointLock = false;
    var refLockValue = 0.0;
    var applyJointHold = false;
    var refHoldValue = 0.0;

    var pumpMode = false;
    var pumpState = 0.0;
    var refractoryTime = 0.0;
    const debounceTime = 0.10 * 2.0 * Math.PI * Math.sqrt(getLength(0) / 9.82) / 2.0;
    var offsetAlpha = 0.0;

    function commonsReset() {
        applyJointLock = false;
        applyJointHold = false;
        pumpMode = false;
    }
    
    function keyDownEvent(e)
    {
        var code = e.keyCode;
        var key = e.key;

        if (key == 'z' || key == 'Z') {
            resetAcrobotState(0.0, 0.0, 0.0, 0.0);
            commonsReset();
        }

        if (key == 'd' || key == 'D') {
            resetAcrobotState(0.01 * (Math.random() - 0.5) - Math.PI / 2, 
                              0.01 * (Math.random() - 0.5) - Math.PI / 2, 
                              0.0, 
                              0.0);
            commonsReset();
        }

        if (key == 'u' || key == 'U') {
            resetAcrobotState(0.01 * (Math.random() - 0.5) + Math.PI / 2, 
                              0.01 * (Math.random() - 0.5) + Math.PI / 2, 
                              0.0, 
                              0.0);
            commonsReset();
        }

        if (key == 'r' || key == 'R') {
            resetAcrobotState(2.0 * Math.PI * Math.random(), 
                              2.0 * Math.PI * Math.random(), 
                              0.0, 
                              0.0);
            commonsReset();
        }

        if (key == 'p' || key == 'P') {
            if (pumpMode) {
                pumpMode = false;
                return;
            }
            //const th1 = -1.0 * Math.PI / 2.0 - 1.0 * Math.random() * Math.PI / 4.0 - Math.PI / 8.0;
            const th1 = -1.0 * Math.PI / 2.0 + 2.0 * (Math.random() - 0.5) * Math.PI / 4.0;
            resetAcrobotState(th1, 
                              th1 - Math.PI / 2.0 + refAngleDelta, 
                              0.0, 
                              0.0);
            refLockValue = getTheta(0) - getTheta(1);
            applyJointLock = true;
            applyJointHold = false;
            offsetAlpha = -1.0 * Math.acos(cosineAlpha(0.0, Math.PI / 2));
            //console.log('offsetAlpha = ' + (offsetAlpha * 180.0 / Math.PI).toFixed(3) + ' deg.');
            pumpState = directionIndicator(offsetAlpha);
            refractoryTime = 0.0;
            pumpMode = true;
        }

        if ((key == 'f' || key == 'F') && !applyJointHold) {
            applyJointLock = !applyJointLock;
            applyTorque(0.0);
            pumpMode = false;
        }

        if ((key == 'h' || key == 'H') && !applyJointLock) {
            applyJointHold = !applyJointHold;
            applyTorque(0.0);
            pumpMode = false;
        }

        if (key == 'b' || key == 'B')
        {
            setFriction(jointWithBrake, 0.50);
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
            if (applyJointLock) {
                refLockValue += refAngleDelta;
            } else if (applyJointHold) {
                refHoldValue += refAngleDelta;
            } else {
                applyTorque(1.0);
            }
        }
        else if (code == 40) // down
        {
            if (applyJointLock) {
                refLockValue -= refAngleDelta;
            } else if (applyJointHold) {
                refHoldValue -= refAngleDelta;
            } else {
                applyTorque(-1.0);
            }
        }
    }

    function keyUpEvent(e)
    {
        var code = e.keyCode;
        var key = e.key;

        if (key == 'b' || key == 'B')
        {
            setFriction(jointWithBrake, 0.0);
        }

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
 
    resetAcrobot();
    resetAcrobotState(0.01 * (Math.random() - 0.5) + Math.PI / 2, 
                      0.01 * (Math.random() - 0.5) + Math.PI / 2, 
                      0.0, 
                      0.0);
    commonsReset();

    console.log(['L1=' + getLength(0), 'L2=' + getLength(1)]);
    console.log([cosineAlpha(0.0, Math.PI / 2), 
                 Math.acos(cosineAlpha(0.0, Math.PI / 2)) * 180 / Math.PI]);
    console.log('debounce time = ' + debounceTime.toFixed(3));

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

        const muA = getFriction(jointWithBrake);
        const muB = getFriction(jointWithFriction);
        ctx.fillText('friction: ' + muA.toFixed(3) + ', ' + muB.toFixed(3) + ' [N m s]', 10.0, 60.0);

        const E = getEnergy();
        ctx.fillText('energy: ' + E.toFixed(3) + ' [J]', 10.0, 80.0);

        if (applyJointLock) {
            var str = 'locked at: ' + refLockValue.toFixed(3) + ' [rad]';
            if (pumpMode) str += ' (pump mode)';
            ctx.fillText(str, 10.0, 100.0);
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
            if (applyJointLock && pumpMode) {
                const presentDirection = directionIndicator(offsetAlpha);
                if (pumpState < 0.0 && presentDirection > 0.0 && refractoryTime > debounceTime) {
                    // "shorten" pendulum
                    refLockValue += 2.0 * refAngleDelta;
                    pumpState = presentDirection;
                    refractoryTime = 0.0;
                }
                if (pumpState > 0.0 && presentDirection < 0.0 && refractoryTime > debounceTime) {
                    // "lengthen" pendulum
                    refLockValue -= 2.0 * refAngleDelta;
                    pumpState = presentDirection;
                    refractoryTime = 0.0;
                }
                refractoryTime += dt;
            }
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
