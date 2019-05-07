% Load parameters from config file. 
 

% Example

%  init:
%    cfg = TextConfig();
%    cfg.declare('mynumber',            3);
%    cfg.declare('mystring', 'teststring');
%    cfg.declare('mynumber2',           3);
%    cfg.loadFile('config.txt');
%
%  file "config.txt" with:
%    mynumber = 2
%    mystring = newstring
%
%  usage:
%    t = cfg.get_num('mynumber')     % -> 2
%    t = cfg.get_str('mystring')     % -> 'newstring'
%    t = cfg.get_num('mynumber2')    % -> 3
%


classdef TextConfig < handle 
   
properties
   % we use this struct as dict (key[str]/value)
   s = struct();
end
   
   
methods
   % constructor 
   function obj = TextConfig()
      obj.s = struct();
   end

   % declare a new parameter
   function declare(self, key, defaultValue)
      assert(ischar(key), 'key has to be a string');
      % always convert value to string
      self.s.(key) = num2str(defaultValue);
   end

   % load parameters from file
   % use declare() first to define parameter names 
   function loadFile(self, filename)
      fprintf('reading config file "%s"\n', filename);
      f = fopen(filename, 'r');
      assert(f >= 0, 'error: unable to open file "%s"', filename);
      while true
         line = fgets(f);
         if (~ischar(line))
            break;
         end
         % remove whitespaces
         line = strtrim(line);
         % skip empty lines
         if (isempty(line))
            continue
         end
         % skip comments
         if (strncmpi(line,'#',1))
            continue
         end
         % split key=value at '='
         t = strsplit(line,'=');
         % t is a 1x2 cell - we have to convert to strings
         key   = strtrim(strjoin(t(:,1)));
         value = strtrim(strjoin(t(:,2)));
         fprintf('loaded %s=%s\n', key, value);
         self.checkField(key);
         self.s.(key) = value;
      end
      fclose(f);
   end

   % is parameter name declared?
   function checkField(self, key)
      assert(isfield(self.s, key), 'config has no member "%s"', key);
   end

   % get string 
   function v = get_str(self, key)
      self.checkField(key);
      v = self.s.(key);
   end

   % get number 
   function v = get_num(self, key)
      v = str2num(self.get_str(key));
      vs = size(v);
      assert(vs(1) == 1, 'config "%s" has size > 1', key)
      assert(vs(2) == 1, 'config "%s" has size > 1', key)
   end

   function v = get_bool(self, key)
      v = str2num(self.get_str(key));
      v = ~(v==0);
   end

   % get vector
   function v = get_vec(self, key)
      v = str2num(self.get_str(key));
      vs = size(v);
      assert(vs(1) == 1, 'config "%s" is not a vector', key)
   end

   % get matrix
   function v = get_mat(self, key)
      value_string = self.get_str(key);
      [v, status] = str2num(value_string);
      assert(status || isempty(value_string), 'config "%s" is not valid: "%s"', key, value_string)
      assert(ismatrix(v), 'config "%s" is not a matrix', key)
   end
   
   % get cell
   function v = get_cell(self, key)
      v = eval(self.get_str(key));
      assert(iscell(v), 'config "%s" is not a cell', key);
   end

   % print all parameters and their values
   function print(self)
      fprintf('configuration:\n');
      keys = fieldnames(self.s);
      for i = 1:numel(keys)
        fprintf('%s\t=\t%s\n', keys{i}, self.s.(keys{i}));
      end
    end
   
end

end
