classdef segTree
    %SEGTREE Class builds and plots a segmentation tree

    % -------------------------------------------------------------
    % Copyright (C) 2021 J.J. Erik Maris
    % Inorganic Chemistry & Catalysis, Utrecht University
    % 
    % This program is free software: you can redistribute it and/or modify
    % it under the terms of the GNU General Public License as published by
    % the Free Software Foundation, either version 3 of the License, or
    % any later version.
    % 
    % This program is distributed in the hope that it will be useful,
    % but WITHOUT ANY WARRANTY; without even the implied warranty of
    % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    % GNU General Public License for more details.
    % 
    % You should have received a copy of the GNU General Public License
    % along with this program.  If not, see <https://www.gnu.org/licenses/>.
    
    properties
        graph
        edgeLabel
        popNodes
    end
    
    methods
        function obj = segTree(nStartNodes)
            % SEGTREE is a constructor of this class
            
            if nargin < 1
                nStartNodes = 1;
            end
            
            obj = obj.init(nStartNodes);
            
        end
        
        function obj = init(obj,nStartNodes)
            
            obj.graph = digraph;
            obj.popNodes = [];
            
            for ii = 1:nStartNodes
                obj.graph = addnode(obj.graph,ii); % add start node
                obj.popNodes = [obj.popNodes; ii];
            end
            
        end
        
        function obj = add(obj,popNo,firstEdgeLabel,secondEdgeLabel,varargin)

            if nargin < 5
                
                popNode = obj.popNodes(popNo); % popnodes has the node number of the each population

                NN = numnodes(obj.graph);
                obj.graph = addedge(obj.graph,[popNode popNode],[NN+1 NN+2]); % make two new nodes
                T = obj.createEdgeTable([popNode; popNode],[NN+1; NN+2],{firstEdgeLabel; secondEdgeLabel});
                obj.edgeLabel = [obj.edgeLabel; T];

                PN = obj.popNodes;
                PN(popNo) = NN+1; % add <=
                obj.popNodes = [PN; NN+2]; % add >
                
            else

                popNode = obj.popNodes(popNo); % popnodes has the node number of the each population

                NN = numnodes(obj.graph);
                NNend = 1:numel(varargin)+2;
                NNend = NNend + NN;
                obj.graph = addedge(obj.graph,repelem(popNode,numel(NNend)),NNend); % make two new nodes
                T = obj.createEdgeTable(repelem(popNode,numel(NNend))',NNend',[{firstEdgeLabel secondEdgeLabel} varargin]');
                obj.edgeLabel = [obj.edgeLabel; T];

                PN = obj.popNodes;
                PN(popNo) = NNend(1); % add first
                obj.popNodes = [PN; NNend(2:end)']; % add other
                
            end
            
        end
        
%         function obj = add(obj,popNo,firstEdgeLabel,secondEdgeLabel)
% 
%             popNode = obj.popNodes(popNo); % popnodes has the node number of the each population
% 
%             NN = numnodes(obj.graph);
%             obj.graph = addedge(obj.graph,[popNode popNode],[NN+1 NN+2]); % make two new nodes
%             obj.edgeLabel = [obj.edgeLabel {firstEdgeLabel secondEdgeLabel}];
% 
%             PN = obj.popNodes;
%             PN(popNo) = NN+1; % add <=
%             obj.popNodes = [PN; NN+2]; % add >
%             
%         end
        
        function plot(obj)

            nodeLabels = cell(numnodes(obj.graph),1);
            nodeLabels(:) = {' '};
            nodeLabels{1} = 'Start';
            for ii = 1:numel(obj.popNodes)
                nodeLabels{obj.popNodes(ii)} = ['Population ' num2str(ii)];
            end
            
            edgeLabels = cell(numedges(obj.graph),1); % preallocate
            idxOut = findedge(obj.graph,obj.edgeLabel.nodeID1,obj.edgeLabel.nodeID2);
            edgeLabels(idxOut) = obj.edgeLabel.nodeName;
            a = plot(obj.graph,'Layout','layered','EdgeLabel',edgeLabels);
            a.NodeLabel = nodeLabels';
        end
        
        function plotEdgeNr(obj)

            nodeLabels = cell(numnodes(obj.graph),1);
            nodeLabels(:) = {' '};
            nodeLabels{1} = 'Start';
            for ii = 1:numel(obj.popNodes)
                nodeLabels{obj.popNodes(ii)} = ['Population ' num2str(ii)];
            end
            
            a = plot(obj.graph,'Layout','layered','EdgeLabel',string(1:numedges(obj.graph)));
            a.NodeLabel = nodeLabels';
        end
        
        function obj = remove(obj,popNo)

            popNode = obj.popNodes(popNo); % popnodes has the node number of the each population

            eid = inedges(obj.graph,popNode); % find edges going to to-be-removed node
            [sOut,tOut] = findedge(obj.graph,eid); % 
            % find idx of desired edge
            bool = (obj.edgeLabel.nodeID1 == sOut & obj.edgeLabel.nodeID2 == tOut);
            obj.edgeLabel(bool,:) = []; % clear labels
            obj.edgeLabel.nodeID1(obj.edgeLabel.nodeID1  > popNode) = obj.edgeLabel.nodeID1(obj.edgeLabel.nodeID1  > popNode) - 1; % restore numebring nodes
            obj.edgeLabel.nodeID2(obj.edgeLabel.nodeID2  > popNode) = obj.edgeLabel.nodeID2(obj.edgeLabel.nodeID2  > popNode) - 1;
            obj.graph = rmnode(obj.graph,popNode); % clear node
            obj.popNodes(obj.popNodes == popNode) = [];
            obj.popNodes(obj.popNodes > popNode) = obj.popNodes(obj.popNodes > popNode) - 1; % update ordering for removed population/node
        end
        
        function obj = merge(obj,popNo1,popNo2)
            
            popNode1 = obj.popNodes(popNo1); % popnodes has the node number of the each population
            popNode2 = obj.popNodes(popNo2); 

            NN = numnodes(obj.graph);
            obj.graph = addedge(obj.graph,[popNode1 popNode2],[NN+1 NN+1]); % make two new nodes
            T = obj.createEdgeTable([popNode1; popNode2],[NN+1; NN+1],{' ';' '}); % don't give these edges a name
            obj.edgeLabel = [obj.edgeLabel; T];

            PN = obj.popNodes;
            PN(popNo1) = NN+1; % keep first population idx
            PN(popNo2) = []; % clear second
            obj.popNodes = PN; % store
            
        end
        
    end
    
    methods(Static)
        function T = createEdgeTable(nodeID1,nodeID2,nodeName)
            % lazy way to construct an edge table
            T = table(nodeID1,nodeID2,nodeName);
        end
    end
end

% function addToDTree(handles)
% S = handles.pushbutton_plotDtree.UserData; % struct
% 
% popNode = S.popNodes(S.popNbr); % popnodes has the node number of the each population
% 
% NN = numnodes(S.graph);
% S.graph = addedge(S.graph,[popNode popNode],[NN+1 NN+2]); % make two new nodes
% 
% if isempty(S.units)
%     unitlabel = [];
% else
%     unitlabel = S.units;
% end
% 
% if isempty(S.edgeLabel)
%     S.edgeLabel = {[S.thName ' <= ' num2str(S.thValue) ' ' unitlabel] [S.thName ' > ' num2str(S.thValue) ' ' unitlabel]};
%     % S.edgeLabel = {[S.thName ' <= ' num2str(S.thValue) ' ' unitlabel] [S.thName ' > ' num2str(S.thValue) ' ' unitlabel]};
% else
%     S.edgeLabel = {S.edgeLabel{:} [S.thName ' <= ' num2str(S.thValue) ' ' unitlabel] [S.thName ' > ' num2str(S.thValue) ' ' unitlabel]};
%     % S.edgeLabel = {S.edgeLabel{:} [S.thName ' <= ' num2str(S.thValue) ' ' unitlabel] [S.thName ' > ' num2str(S.thValue) ' ' unitlabel]};
% end
% 
% PN = S.popNodes;
% PN(S.popNbr) = NN+1; % add <=
% S.popNodes = [PN; NN+2]; % add >
% 
% handles.pushbutton_plotDtree.UserData = S;
% 
% function deleteFromDTree(handles,popNo)
% S = handles.pushbutton_plotDtree.UserData; % struct
% 
% popNode = S.popNodes(popNo); % popnodes has the node number of the each population
% 
% eid = inedges(S.graph,popNode); % find edges going to to-be-removed node
% S.edgeLabel(eid) = []; % clear labels
% S.graph = rmnode(S.graph,popNode); % clear node
% S.popNodes(popNode) = [];
% 
% handles.pushbutton_plotDtree.UserData = S;
% 
% function plotDTree(handles)
% 
% S = handles.pushbutton_plotDtree.UserData; % struct
% nodeLabels = cell(numnodes(S.graph),1);
% nodeLabels(:) = {''};
% nodeLabels{1} = 'Start';
% for ii = 1:numel(S.popNodes)
%     nodeLabels{S.popNodes(ii)} = ['Population ' num2str(ii)];
% end
% figure
% a = plot(S.graph,'Layout','layered','EdgeLabel',S.edgeLabel);
% a.NodeLabel = nodeLabels';
