function [ phiPrime ] = backwardsDifference( delT, phi, phiMinusOne, phiMinusTwo, phiMinusThree, phiMinusFour, phiMinusFive, phiMinusSix, order )
%backwardsDifference: Calculates the derivative based on
%previously known value. This is either first through sixth order accurate.
    switch order
        case 1
            phiPrime = (phi - phiMinusOne) ./ delT;
        case 2
            phiPrime = (3.*phi - 4.*phiMinusOne + phiMinusTwo) ./ (2*delT);
        case 3
            phiPrime = ((11/16)*phi - 3*phiMinusOne + (3/2)*phiMinusTwo - (1/3)*phiMinusThree) / delT;
        case 4
            phiPrime = ((25/12)*phi - 4*phiMinusOne + 3*phiMinusTwo - (4/3)*phiMinusThree + (1/4)*phiMinusFour) / delT;
        case 5
            phiPrime = ((137/60)*phi - 5*phiMinusOne + 5*phiMinusTwo - (10/3)*phiMinusThree + (5/4)*phiMinusFour - (1/5)*phiMinusFive) / delT;
        case 6
            phiPrime = ((49/20)*phi - 6*phiMinusOne + (15/2)*phiMinusTwo - (20/3)*phiMinusThree + (15/4)*phiMinusFour - (6/5)*phiMinusFive + (1/6)*phiMinusSix) / delT;
        otherwise
            fprintf('ERROR! ''%i'' is an invalid choice for order. Valid choices are ''1,'' ''2,'' ''3,'' ''4,'' ''5,'' ''6.''',order)
    end
end

