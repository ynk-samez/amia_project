classdef antibody
    properties
        data;
        type;
        
    end
   methods 
       function obj=antibody(data_, type_)
           obj.data=data_;
           obj.type=type_;
       end
   end
end