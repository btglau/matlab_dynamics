if ~strcmp(getenv('PBS_JOBID'),''); 
%If running on cluster

    delete(gcp('nocreate'))
    
        fprintf('\nAttempting to remove cluster storage directory...');
    [status, message, messageid] = rmdir(local_clust_data,'s');
    if status
        fprintf('success.\n');
    else
        fprintf(['\n\nError - could not remove %s for some reason.\n'...
            'Looks like it needs to be manually removed.\n Error: %s, %s\n'],local_clust_data,message,messageid);
    end
end