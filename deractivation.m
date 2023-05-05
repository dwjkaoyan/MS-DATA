function [derneuron_ot] = deractivation(neuron_state,Neuron_type,slope)
%Neuron_type: 1-sigmoid; 2- tanh
X=slope*neuron_state;
if Neuron_type==1    
%     derneuron_ot=exp(-X)./(1+exp(-X)).^2;
    derneuron_ot=logsig(slope*X).*(1-logsig(slope*X));
else derneuron_ot=1-(tanh(slope*X)).^2;
end

    