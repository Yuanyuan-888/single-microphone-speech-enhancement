
function [SNR,Pss] = SpeechPSDML(Pyy,Pnn)
SNR = Pyy ./Pnn - 1; 
Pss = Pyy - Pnn;
end
