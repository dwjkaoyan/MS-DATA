function [neuron_ot] = activation(neuron_state,Neuron_type,slope)
%Neuron_type: 1-sigmoid; 2- tanh
if Neuron_type==1
    neuron_ot = logsig(slope*neuron_state);
else if Neuron_type==2 
        neuron_ot = tanh(slope*neuron_state);
    else neuron_ot = neuron_state;
    end
end
end