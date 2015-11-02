NodeNums = 100; % the num of node  
AreaR = 100 ;   % the area of simulate 
NodeTranR=10;   % the transit Radius 
Elec=50 * 10^(-9); % 
Eamp=100*10^(-12);  
Bx=50; % The Postion of Baseation 
By=175; 
MaxInteral =500; % the leach simulate time 
Pch=0.05; % the desired percentage of cluster heads  
InitEn=0.5; % the init energy of all node 
Tr=30;  
TDMA=1; 
Kbit=2000; % the bits of a node transmiting a packet every time 
BandWitch = 1*10.^(6); % Channel Bandwitch 
TOS_LOCAL_ADDRESS = 0;  
for i=1:(MaxInteral) 
    AliveNode(i)=NodeNums; 
    AmountData(i)=0; 
end 
sym alldata; 
alldata=0; 
LAECH = zeros(1,MaxInteral); 
LAENO = zeros(1,MaxInteral);  
for i=1:1:NodeNums  
    EnNode(i)=InitEn; % the init energy of all node 
    StateNode(i)=1;    % the State of all node 1: alive 0:dead 
    ClusterHeads(i)=0; % the Set of Cluster Head ,1: cluster head 0 :node 
     
    Rounds=0; % the round 
end 
 
 
Threshold=0;    % the threshold of node becoming a cluster-head 
 
    Node.x=AreaR*rand(1,NodeNums); % the position of node  
    Node.y=AreaR*rand(1,NodeNums); 
    Node.c=zeros(1,NodeNums); 
    Node.d=zeros(1,NodeNums); 
    Node.l=zeros(1,NodeNums); 
    Node.csize=zeros(1,NodeNums); 
    Node.initclEn=zeros(1,NodeNums); 
%    for i=1:NodeNums 
%     Node.c(i)=0;                   % the Cluster head of node  
%     Node.d(i)=0;                   % the distance between cluster head and node 
%     Node.l(i)=Kbit;                % the length of node i transmit packet 
%     Node.csize(i)=0; 
%    end 
 
for Rounds = 1:MaxInteral  
   % the Setup phase of cluster 
    Node.csize=Node.csize-Node.csize; 
    Node.d=Node.d-Node.d; 
    Node.c=Node.c-Node.c;  
    for i =1:NodeNums 
      Threshold=Pch/(1-Pch*(mod(Rounds-1,1/Pch))); 
%         Threshold=(Pch/(1-Pch*(mod(Rounds-1,1/Pch))))*(Node.initclEn/InitEn); 
% new=((countch)/(1/Pch))*(1-(Node.initclEn/InitEn));
%           Threshold=(Pch/(1-Pch*(mod(Rounds-1,1/Pch))))*((Node.initclEn/InitEn)+new); 
      if StateNode(i)==1         % if node is alive 
          if ClusterHeads(i) ==1  
              ClusterHeads(i)=0; 
          elseif rand(1,1)<Threshold 
           
          ClusterHeads(i)=1; 
          Node.c(i)=TOS_LOCAL_ADDRESS; 
           Node.initclEn(i)=EnNode(i); 
         else   ClusetHeads(i)=0; 
                Node.initclEn(i)=EnNode(i); 
         end 
      end 
    end 
   if sum(ClusterHeads)==0 
       continue; 
   end 
     
    EntranPCH = Elec * Kbit+ Eamp*Kbit*((Tr.^2+Tr.^2)); % The expended engergy by new Cluster head advertising that it is new cluster head 
    for i=1:NodeNums      
        if ClusterHeads(i) ==1  
          
          if EnNode(i) >= EntranPCH 
                  EnNode(i) = EnNode(i) - EntranPCH ; 
              else  
                  StateNode(i)=0; 
           end    
       end  
    end 
     
    for i=1:NodeNums 
      if StateNode(i)==1         % if node is alive 
        if ClusterHeads(i) ~=1   % the node is not cluster head 
          for j=1:NodeNums        
             if ClusterHeads(j) ==1    
               dist = ((Node.x(i)-Node.x(j)).^2)+((Node.y(i)-Node.y(j)).^2); % the distance.^2 
                EnRecP = Elec * Kbit ; 
                  if EnNode(i) >= EnRecP            % the energy reciving a boardcast packet can expend     
                    EnNode(i) = EnNode(i) - EnRecP ; 
                  else  
                    StateNode(i)=0; 
                  end     
                  if Node.d(i) ==0                  % choose the cluster head 
                   Node.d(i)=dist ; 
                   Node.c(i)=j;  
                  else 
                   if Node.d(i) > dist 
                       Node.d(i)=dist ; 
                       Node.c(i)=j;  
                   end 
                  end 
             end 
         %%%%% end of choosing the cluster head ,Node.c(i) save the id of 
         %%%%% cluster head 
         end 
          
           if StateNode(i)==1 
            Node.csize(Node.c(i))= Node.csize(Node.c(i))+1; 
           end 
        else      % the node is cluster head 
         Node.d(i)=((Node.x(i)-Bx).^2)+((Node.y(i)-By).^2) ; 
         Node.c(i)=TOS_LOCAL_ADDRESS;  
        end 
     end 
    end 
     
   % painting the node and the cluster head 
   % for i=1:NodeNums 
   %     if ClusterHeads(i)==1 
   %         plot(Node.x(i),Node.y(i),'rs'); 
   %         hold on; 
   %     else plot(Node.x(i),Node.y(i),'k*');    
   %           hold on; 
   %     end 
   % end 
    
    % the TDMA Phase 
    alldata=0; 
      for i=1:NodeNums 
        if StateNode(i)==1  
         if ClusterHeads(i)==1 
              
             TolLengthPacket = Kbit.*Node.csize(i); 
             alldata=alldata+TolLengthPacket; 
             EntranPCH = Elec * TolLengthPacket+ Eamp*TolLengthPacket*(Node.d(i)); 
             EntranPCH.*TDMA; 
              if EnNode(i) >= EntranPCH 
                   EnNode(i) = EnNode(i) - EntranPCH ; 
               else  
                   StateNode(i)=0; 
              end  
          else 
               EntranP = Elec *Node.l(i)+ Eamp*Node.l(i)*(Node.d(i));  
               EntranP=EntranP.*TDMA; 
               if EnNode(i) >= EntranP  
                   EnNode(i) = EnNode(i)-EntranP; 
               else 
                   StateNode(i)=0;    % the node dead  
               end 
                EnRecP = Elec * Node.l(i) ; 
                EnRecP=EnRecP.*TDMA; 
               if EnNode(Node.c(i)) >= EnRecP 
                   EnNode(Node.c(i)) = EnNode(i) - EnRecP ; 
               else  
                   StateNode(Node.c(i))=0; 
               end  
          end 
        end  
      end 
   if Rounds==1 
       AmountData(Rounds)=alldata; 
   else 
      AmountData(Rounds)=alldata+AmountData(Rounds-1); 
   end 
 
    for i=1:NodeNums  
         if StateNode(i)==0 
             AliveNode(Rounds)= AliveNode(Rounds)-1; 
         end 
     end       
   % the TDMA Phase       
%    for RNum=1:TDMA 
%     for i=1:NodeNums 
%       if ClusterHeads(i)==1 
%          % EntranPCH = Elec * Node.l(i)+ Eamp*Node.l(j)*(Eamp*Node.l(j).^2);  
%          TolLengthPacket = 0; 
%           for j=1:NodeNums 
%            if Node.c(j) ==i 
%               TolLengthPacket =TolLengthPacket+ Node.l(j); 
%               EntranP = Elec * Node.l(j)+ Eamp*Node.l(j)*(Node.d(j)); % The require energy of node transmitting a packet  
%               EnRecP = Elec * Node.l(j) ; % The require energy of recving a packet  
%               if EnNode(j) >= EntranP  
%                   EnNode(j) = EnNode(j)-EntranP; 
%               else 
%                   StateNode(j)=0;    % the node dead  
%               end 
%               if EnNode(i) >= EnRecP 
%                   EnNode(i) = EnNode(i) - EnRecP ; 
%               else  
%                   StateNode(i)=0; 
%               end     
%            end 
%           end 
%           EntranPCH = Elec * TolLengthPacket+ Eamp*TolLengthPacket*(Node.d(j)); 
%           if EnNode(i) >= EntranPCH 
%                   EnNode(i) = EnNode(i) - EntranPCH ; 
%               else  
%                   StateNode(i)=0; 
%               end     
%       end 
%     end 
% end 
syms sumch sumno countch countno ; 
sumch=0;  
sumno=0; 
countch=0; 
countno=0; 
    for i=1:NodeNums  
        if Node.initclEn(i)>0 
           if ClusterHeads(i)==1 
              sumch=sumch+Node.initclEn(i); 
              countch=countch+1; 
           else 
              sumno=sumno+ Node.initclEn(i); 
              countno=countno+1; 
           end 
         end 
        %if StateNode(i)==0 
        %    AliveNode(Rounds)= AliveNode(Rounds)-1; 
        %    if AliveNode(Rounds) == 0 
        %        break; 
        %    end 
        %end 
    end 
    LAECH(Rounds) = sumch/countch; 
    LAENO(Rounds) = sumno/countno;  
    
   % Rounds = Rounds+1 
end 
xtime= 1:1:MaxInteral; 
plot(xtime,AliveNode);   


% 
%clear ;  
   % plot(Node.x,Node.y); 
% plot(Node.x,Node.y,'--rs','LineWidth',2,... 
%              'MarkerEdgeColor','k',... 
%              'MarkerFaceColor','g',... 
%              'MarkerSize',10); 
% plot(Node.x,Node.y,'k*');