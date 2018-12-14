classdef PottsModel < handle
    %POTTSMODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        state
        HSV
        m
        n
        ax_handle
    end
    
    methods
        function obj = PottsModel(ax_handle,m,n,N)
            %POTTSMODEL Construct an instance of this class
            %   Detailed explanation goes here
            assert(m > 0);
            assert(n > 0);
            obj.m = m;
            obj.n = n;
            obj.HSV = ones(m,n,3);
            obj.HSV(:,:,2) = 0.95;
            obj.HSV(:,:,3) = 0.75;
            obj.randomize(N);
            assert(isa(ax_handle,'handle'));
            obj.ax_handle = ax_handle;
            obj.plot();
        end
        
        function randomize(obj,N)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            thetas = linspace(-pi,pi,N+1);
            obj.state = randi(N,obj.m,obj.n);
            for i = 1:N
                obj.state(obj.state == i) = thetas(i);
            end
        end
        
        function h = plot(obj)
            obj.HSV(:,:,1) = obj.state/(2*pi) + 0.5;
            h = imshow(hsv2rgb(obj.HSV),'Parent',obj.ax_handle);
        end
        
        function MC_step(obj,N,J,h,h_theta,alpha)
            old_state = obj.state;
            theta2s = linspace(-pi,pi,N+1);
            for i = 1:obj.m
                for j = 1:obj.n
                    k = ceil(rand()*N);
                    theta2 = theta2s(k);
                    if obj.state(i,j) ~= theta2
                        dH = -obj.Energy(i,j,theta2,J,h,h_theta) + ...
                                obj.Energy(i,j,old_state(i,j),J,h,h_theta);
                        if dH < 0
                            obj.state(i,j) = theta2;    
                        elseif exp(-dH*alpha) > rand()
                            obj.state(i,j) = theta2;
                        end
                    end
                end
            end
        end
        
        function kMC_step(obj,N,J,h,h_theta,alpha)
            energy = zeros(obj.m,obj.n,N);
            theta2s = linspace(-pi,pi,N+1);
            for i = 1:obj.m
                for j = 1:obj.n
                    for k = 1:length(theta2s)-1
                        theta2 = theta2s(k);
                        if obj.state(i,j) == theta2s(k)
                            energy(i,j,k) = 0;
                        else
                            energy(i,j,k) = ...
                                -obj.Energy(i,j,theta2,J,h,h_theta) + ...
                                obj.Energy(i,j,obj.state(i,j),J,h,h_theta);
                        end
                    end
                end
            end
            omega = exp(-alpha*energy);
            omega_tot = sum(sum(sum(omega)));
            r1 = omega_tot*rand();
            running_tot = 0;
            b = false;
            for i = 1:obj.m
                for j = 1:obj.n
                    for k = 1:length(theta2s)-1
                        running_tot = running_tot + omega(i,j,k);
                        if running_tot >= r1
                            obj.state(i,j) = theta2s(k);
                            b = true;
                        end
                        if b
                            break;
                        end
                    end
                    if b
                        break;
                    end
                end
                if b
                    break;
                end
            end
        end
        
        function H = Energy(obj,i,j,theta,J,h,h_theta)
            NNij = obj.neighbor_indecies(i,j);
            H = 0;
            for i = 1:length(NNij(:,1))
                H = H + J*cos(theta - obj.state(NNij(i,1),NNij(i,2))) + ...
                    h*cos(theta - h_theta);
            end
            return;
        end
        
        function NNij = neighbor_indecies(obj,i,j)
            if i == 1
                up = [obj.m,j];
            else
                up = [i-1,j];
            end
            if i == obj.m
                down = [1,j];
            else
                down = [i+1,j];
            end
            if j == 1
                left = [i,obj.n];
            else
                left = [i,j-1];
            end
            if j == obj.n
                right = [i,1];
            else
                right = [i,j+1];
            end
            NNij = [up;down;left;right];
            return;
        end
        
    end
end

