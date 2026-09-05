% Rebuild every P-file in the project from its source
%
% Author : Gustav Strijkers
% Date   : 2026-09-05

function pcodeAll(varargin)

    % Purpose:
    %   Rebuild the project's P-files from their .m sources, all with the release
    %   flag retro.pcodeTarget gives, so no file is left in a format the rest of
    %   the tree does not share.
    %
    % How it works:
    %   Every .m that already has a .p beside it is rebuilt; a source with no P-file
    %   is left alone, so what is p-coded stays a property of the tree rather than a
    %   list kept here. Sources under tests/ are skipped, since the suite runs from
    %   source.
    %
    %   pcode writes into a +package/ layout under the working directory and cannot
    %   be trusted to work in place on a class inside a package, so each file is
    %   built in a temporary folder and the result copied over its .p. The temporary
    %   folder is removed afterwards whatever happens.
    %
    %   An optional filter builds only the files whose path contains it, for working
    %   on one folder.
    %
    % Inputs:
    %   filter - char row, optional, a substring of the paths to rebuild
    %
    % Output:
    %   none - the P-files are written in place, and a count is printed
    %
    % Raises:
    %   Retrospective:pcodeAll - a source failed to compile

    root = fileparts(mfilename('fullpath'));
    flag = retro.pcodeTarget();

    filter = '';
    if nargin > 0
        filter = varargin{1};
    end

    sources = dir(fullfile(root, '**', '*.m'));
    tmp = fullfile(tempdir, ['pcodeAll_' char(java.util.UUID.randomUUID)]);
    mkdir(tmp);
    cleanup = onCleanup(@() rmdir(tmp, 's'));

    built = 0;
    skipped = 0;

    for k = 1:numel(sources)

        source = fullfile(sources(k).folder, sources(k).name);
        target = replace(source, '.m', '.p');

        if contains(source, [filesep 'tests' filesep]) || ~isfile(target)
            skipped = skipped + 1;
            continue
        end

        if ~isempty(filter) && ~contains(source, filter)
            skipped = skipped + 1;
            continue
        end

        here = cd(tmp);
        try
            pcode(source, flag);
        catch ME
            cd(here);
            error('Retrospective:pcodeAll', 'pcode failed on %s: %s', source, ME.message);
        end
        cd(here);

        [~, name] = fileparts(source);
        made = dir(fullfile(tmp, '**', [name '.p']));

        if isempty(made)
            error('Retrospective:pcodeAll', 'pcode wrote nothing for %s', source);
        end

        copyfile(fullfile(made(1).folder, made(1).name), target);
        delete(fullfile(made(1).folder, made(1).name));
        built = built + 1;

    end

    fprintf('pcodeAll: %d P-files rebuilt with %s, %d sources skipped\n', built, flag, skipped);

end % pcodeAll
