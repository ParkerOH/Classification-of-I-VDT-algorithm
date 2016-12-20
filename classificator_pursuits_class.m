% I-DT model classificator
classdef classificator_pursuits_class <  eye_tracker_raw_data_reader_class & ...             % Reader from eye tracker data
                                         eye_records_class & ...                             % Basic class for placing eye tracker data
                                         eye_tracker_raw_data_converter_ETU_degree & ...     % Convertor between ETU and degress in data
                                         eye_tracker_raw_data_filter_class & ...             % Eye tracker data filtering by range of degrees
                                         classificator_merge_class & ...                     % Creates sequences of eye movements
                                         classificator_saccade_amplitude_filter_class & ...  % Filtered saccades based theire amplitude
                                         classificator_datafile_output_class & ...           % Output sequences to the files
                                         classificator_get_percentage_class & ...            % Calculate percentage of movements of every type
                                         classificator_enumerations_class & ...              % Basic enumerations definitions
                                         classificator_time_rate_class & ...                 % Time step and sample rate definitions
                                    handle
    % This is skeleton class for user classification
  
    properties(Hidden)
  %      saccade_detection_threshold;
  %      dispersion_threshold;
    end

    methods

% Classification function
        function classify(obj)
                  
            if( obj.debug_mode ~= 0)
                fprintf(strcat('Begin data classification with user classifier in :',datestr(now),'\n'));
            end
            obj.calculate_delta_t();
            saccade_detection_threshold = 70;
            dispersion_threshold = .8;
            x_velocity_degree = zeros( length(obj.eye_records),1 );
            y_velocity_degree = zeros( length(obj.eye_records),1 );
% Calculate absolute degree velocity of our records
            x_velocity_degree( 2:end ) =(   obj.eye_records( 2:end,obj.X_COORD ) - ...
                                            obj.eye_records( 1:end-1,obj.X_COORD ) ) / obj.delta_t_sec;
            y_velocity_degree( 2:end ) =(   obj.eye_records( 2:end,obj.Y_COORD ) - ...
                                            obj.eye_records( 1:end-1,obj.Y_COORD ) ) / obj.delta_t_sec;
% First point is a special case
            x_velocity_degree(1) = 0;
            y_velocity_degree(1) = 0;
            obj.eye_records(:,obj.VELOCITY) = sqrt( x_velocity_degree.^2 + y_velocity_degree.^2 );
% First point is a special case
            obj.eye_records(1,obj.MOV_TYPE ) = obj.NOISE_TYPE;
            obj.eye_records(1,obj.VELOCITY) = 0;
% Now we mark fixations            
            obj.eye_records( ((abs(obj.eye_records(:,obj.VELOCITY)) < saccade_detection_threshold) & (obj.eye_records(:,obj.MOV_TYPE) ~= obj.NOISE_TYPE)),obj.MOV_TYPE ) = obj.FIXATION_TYPE;
% Now we mark saccades
            obj.eye_records( ((abs(obj.eye_records(:,obj.VELOCITY)) >= saccade_detection_threshold) & (obj.eye_records(:,obj.MOV_TYPE) ~= obj.NOISE_TYPE)),obj.MOV_TYPE ) = obj.SACCADE_TYPE;
% Now we mark every invalid point as noise
            obj.eye_records( (obj.eye_records(:,obj.VALIDITY) == obj.DATA_INVALID),obj.MOV_TYPE ) = obj.NOISE_TYPE;
% And we mark every noise points as invalid
            obj.eye_records( (obj.eye_records(:,obj.MOV_TYPE) == obj.NOISE_TYPE),obj.VALIDITY ) = obj.DATA_INVALID;
% Now we have to expand saccade mark to previous point
            tmp_type = obj.eye_records(2:end,obj.MOV_TYPE);
            tmp_type(length(tmp_type)+1) = NaN;
            obj.eye_records( ((obj.eye_records(:,obj.VALIDITY) == obj.DATA_VALID) & ...
                              (tmp_type(:) == obj.SACCADE_TYPE)),obj.MOV_TYPE) = obj.SACCADE_TYPE;
            obj.eye_records( 1 , obj.MOV_TYPE ) = obj.eye_records( 2 , obj.MOV_TYPE );
%First we read in the record array, and its parts. 
%We declare counter as a counter, a flag to end the loop, and position as
%   the location in the array.
%
          
            size =0;
            for i =1: length(obj.eye_records)
                if(obj.eye_records(i,obj.MOV_TYPE) ~= obj.SACCADE_TYPE && obj.eye_records(i,obj.NOISE_TYPE))
                   size = size + 1;
                   record(size).ID = i;
                   record(size).X_COORD = obj.eye_records(i,obj.X_COORD);
                   record(size).T_COORD = obj.eye_records(i,obj.T_COORD);
                   record(size).Y_COORD = obj.eye_records(i,obj.Y_COORD);
                end
            end
            counter = 1;
            flag = false;
            position = 1;
            while(counter < length(record) && flag == false)
                tempw(1) = record(counter);
                if(counter ~= length(record))
                   tempw(2) = record(counter + 1);
                   position = 2;
                else        
                   flag = true;
                   obj.eye_records(record(length(record)).ID ,obj.MOV_TYPE ) = obj.PURSUIT_TYPE;
                   position = 1;
                end               
                counter = counter + 1;
%Here from the temporary window array variable, tempw, we read in the x and y variables
%Then we compare for min and max
                x_Max = tempw(1).X_COORD;
                x_Min = tempw(1).X_COORD;
                y_Max = tempw(1).Y_COORD;
                y_Min = tempw(1).Y_COORD;
                for i=1:length(tempw)
                    if (tempw(i).X_COORD > x_Max)
                        x_Max = tempw(i).X_COORD;
                    end
                    if(tempw(i).X_COORD < x_Min)
                        x_Min = tempw(i).X_COORD;
                    end
                    if(tempw(i).Y_COORD > y_Max)
                        y_Max = tempw(i).Y_COORD;
                    end
                    if (tempw(i).Y_COORD < y_Min)
                        y_Min = tempw(i).Y_COORD;
                    end
                end
                dispersion = (abs(x_Max - x_Min) + abs(y_Max - y_Min));    
%Now that our dispersion is made we can make our window and scan the points
                if(dispersion < dispersion_threshold)
                    while(dispersion < dispersion_threshold && counter < length(record))
                        tempw(position+1) = record(counter+1);
                        for i=1:length(tempw)
                            if (tempw(i).X_COORD > x_Max)
                                x_Max = tempw(i).X_COORD;
                            end
                            if(tempw(i).X_COORD < x_Min)
                               x_Min = tempw(i).X_COORD;
                            end
                            if(tempw(i).Y_COORD > y_Max)
                               y_Max = tempw(i).Y_COORD;
                            end
                            if (tempw(i).Y_COORD < y_Min)
                                y_Min = tempw(i).Y_COORD;
                            end
                        end
                        dispersion = (abs(x_Max - x_Min) + abs(y_Max - y_Min));
                        position = position + 1;
                        counter  = counter + 1;
                    end
                else
                    obj.eye_records(tempw(1).ID ,obj.MOV_TYPE ) = obj.PURSUIT_TYPE;
                end
            end
            if(obj.debug_mode ~= 0)
                   fprintf(strcat('Complete data classification with user classifier in :',datestr(now),'\n'));
       
            end
        end
        %function set.saccade_detection_threshold(obj,value)
        %         obj.saccade_detection_threshold = value;
        %end
    end
end
